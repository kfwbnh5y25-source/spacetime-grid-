import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { execFileSync } from 'node:child_process';

const rootDir = path.resolve(import.meta.dirname, '..');
const indexHtml = fs.readFileSync(path.join(rootDir, 'index.html'), 'utf8');
const feedbackHtml = fs.readFileSync(path.join(rootDir, 'feedback-lab', 'index.html'), 'utf8');
const swJs = fs.readFileSync(path.join(rootDir, 'sw.js'), 'utf8');

function extractInlineScript(html, label) {
  const matches = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)];
  assert.equal(matches.length, 1, `Expected 1 inline script in ${label}, found ${matches.length}`);
  return matches[0][1];
}

const siteScript = extractInlineScript(indexHtml, 'index.html');
const feedbackScript = extractInlineScript(feedbackHtml, 'feedback-lab/index.html');

function loadPhysicsContext() {
  const noop = () => {};
  const makeCtx = () => new Proxy({}, { get: () => noop });
  const makeEl = () => ({
    style: {},
    classList: { add: noop, remove: noop },
    setAttribute: noop,
    addEventListener: noop,
    getContext: (type) => (type === '2d' ? makeCtx() : null),
    querySelector: () => makeEl(),
    textContent: ''
  });
  const sandbox = {
    console,
    Math,
    Number,
    Array,
    Map,
    Float32Array,
    Float64Array,
    Uint32Array,
    URL,
    URLSearchParams,
    Blob,
    Image: class { constructor() {} },
    PerformanceObserver: undefined,
    innerWidth: 1200,
    innerHeight: 800,
    devicePixelRatio: 2,
    screen: { width: 1200, height: 800 },
    location: { search: '', hash: '', pathname: '/', protocol: 'https:', hostname: 'localhost', href: 'https://localhost/' },
    navigator: { userAgent: 'node-test' },
    localStorage: { getItem: () => null, setItem: noop },
    history: { replaceState: noop },
    matchMedia: () => ({ matches: false, addEventListener: noop }),
    performance: { now: () => 0, getEntriesByName: () => [] },
    requestAnimationFrame: () => 0,
    cancelAnimationFrame: noop,
    setTimeout: () => 0,
    clearTimeout: noop,
    addEventListener: noop,
    document: {
      hidden: true,
      activeElement: null,
      getElementById: () => makeEl(),
      createElement: () => makeEl(),
      querySelectorAll: () => [],
      addEventListener: noop,
      body: { appendChild: noop }
    },
    window: null
  };
  sandbox.window = sandbox;

  const instrumented = siteScript.replace(
    /\}\)\(\);\s*$/,
    `  globalThis.__physics = {
    DISPLAY_SCHWARZSCHILD_SCALE,
    STATIC_WELL_PULL_PX,
    STATIC_WELL_SOFTENING_PX,
    CAUSAL_FRONT_WIDTH_PX,
    SLINGSHOT_SPEED_SCALE,
    N_BODY_MU_PER_MASS,
    MAX_GPU_POINTERS,
    staticSchwarzschildDisplacement,
    deformSpacetimePoint,
    runGeodesicValidation,
    geodesicValidation,
    newMassValue,
    gravityG,
    gridCameraState,
    worldPoint,
    screenPoint,
    persistentMasses,
    createMass,
    massAccelerationAt,
    advanceMasses,
    classifyLaunchOrbit,
    predictLaunchArc,
    resetView,
    get cam() { return { viewX, viewY, viewScale, viewRotation, camComX, camComY, width, height }; },
    setCam(patch) {
      if ('viewX' in patch) viewX = patch.viewX;
      if ('viewY' in patch) viewY = patch.viewY;
      if ('viewScale' in patch) viewScale = patch.viewScale;
      if ('viewRotation' in patch) viewRotation = patch.viewRotation;
      if ('camComX' in patch) camComX = patch.camComX;
      if ('camComY' in patch) camComY = patch.camComY;
    }
  };
})();`
  );

  const context = vm.createContext(sandbox);
  vm.runInContext(instrumented, context);
  return context.__physics;
}

const physics = loadPhysicsContext();

test('A. Static curvature: prominent, symmetric, smooth, and time-independent after front passes', () => {
  const mass = { x: 600, y: 400, mass: physics.newMassValue, createdAt: 0, strength: 1 };
  // Once the causal front at c=180 px/s has passed all sample points (t >= 5000 ms),
  // curvature must be identical across all later timestamps.
  const timestamps = [5000, 10000, 30000, 60000];
  const samplePoints = [
    { x: 655, y: 400 },
    { x: 545, y: 400 },
    { x: 600, y: 455 },
    { x: 600, y: 345 },
    { x: 720, y: 400 },
    { x: 900, y: 400 }
  ];

  const baseline = samplePoints.map((pt) => physics.staticSchwarzschildDisplacement(pt.x, pt.y, timestamps[0], [mass]));
  for (const t of timestamps) {
    samplePoints.forEach((pt, idx) => {
      const d = physics.staticSchwarzschildDisplacement(pt.x, pt.y, t, [mass]);
      assert.equal(d.x, baseline[idx].x, `Static dx changed at t=${t}`);
      assert.equal(d.y, baseline[idx].y, `Static dy changed at t=${t}`);
    });
  }

  // Radial symmetry around the stationary mass (55 px in +x, -x, +y, -y)
  const dPx = physics.staticSchwarzschildDisplacement(655, 400, 10000, [mass]);
  const dNx = physics.staticSchwarzschildDisplacement(545, 400, 10000, [mass]);
  const dPy = physics.staticSchwarzschildDisplacement(600, 455, 10000, [mass]);
  const dNy = physics.staticSchwarzschildDisplacement(600, 345, 10000, [mass]);
  assert.ok(dPx.x < 0 && Math.abs(dPx.y) < 1e-12, 'Pull on +x axis points toward mass (-x)');
  assert.ok(Math.abs(dPx.x + dNx.x) < 1e-12, 'Symmetric between +x and -x');
  assert.ok(Math.abs(dPx.x - dPy.y) < 1e-12, 'Symmetric between x and y axes');
  assert.ok(Math.abs(dPy.y + dNy.y) < 1e-12, 'Symmetric between +y and -y');

  // Prominent visibility (~20-25 px peak pull near mass, decaying smoothly with distance)
  const pull55 = Math.hypot(dPx.x, dPx.y);
  const pull120 = Math.hypot(...Object.values(physics.staticSchwarzschildDisplacement(720, 400, 10000, [mass])));
  const pull300 = Math.hypot(...Object.values(physics.staticSchwarzschildDisplacement(900, 400, 10000, [mass])));
  assert.ok(pull55 >= 18 && pull55 <= 28, `Expected prominent 18-28 px pull at 55 px, got ${pull55.toFixed(2)} px`);
  assert.ok(pull120 > 12 && pull120 < pull55, `Expected smooth decay at 120 px, got ${pull120.toFixed(2)} px`);
  assert.ok(pull300 < pull120 && pull300 > 5, `Expected continued 1/r decay at 300 px, got ${pull300.toFixed(2)} px`);
});

test('B. Causal curvature propagation: flat ahead of light cone, smooth turn-on at c, persistent well behind', () => {
  const mass = { x: 600, y: 400, mass: physics.newMassValue, createdAt: 0, strength: 1 };

  // At r = 360 px, causal front center arrives at t = 360/180 = 2.0 s (leading edge at 2.0 - 28/180 = 1.844 s).
  // Ahead of the leading edge (t = 1500 ms), space must be completely flat.
  const preFront = physics.staticSchwarzschildDisplacement(960, 400, 1500, [mass]);
  assert.equal(preFront.x, 0, 'Curvature must be zero ahead of the causal front');
  assert.equal(preFront.y, 0, 'Curvature must be zero ahead of the causal front');

  // At t = 2000 ms (front center right at r = 360 px), activation = 0.5
  const midFront = physics.staticSchwarzschildDisplacement(960, 400, 2000, [mass]);
  const fullStatic = physics.staticSchwarzschildDisplacement(960, 400, 10000, [mass]);
  assert.ok(Math.abs(midFront.x - 0.5 * fullStatic.x) < 1e-6, 'At front center, curvature should be 50% activated');

  // Behind the front (t = 2500 ms and t = 30000 ms), curvature stays permanently at full static value (no wave oscillation)
  const postFront = physics.staticSchwarzschildDisplacement(960, 400, 2500, [mass]);
  const lateStatic = physics.deformSpacetimePoint(960, 400, 30000, [mass]);
  assert.ok(Math.abs(postFront.x - fullStatic.x) < 1e-6, 'Full curvature reached behind causal front');
  assert.ok(Math.abs((lateStatic.x - 960) - fullStatic.x) < 1e-6, 'Curvature persists unchanged at late times');
});

test('C. Mass scaling: monotonic curvature increase across 3 mass strengths and no line folding', () => {
  const masses = [0.001, 0.002, 0.006];
  const pulls = masses.map((m) => {
    const src = { x: 600, y: 400, mass: m, createdAt: 0, strength: 1 };
    const d = physics.staticSchwarzschildDisplacement(655, 400, 10000, [src]);
    return Math.hypot(d.x, d.y);
  });
  assert.ok(pulls[0] < pulls[1] && pulls[1] < pulls[2], `Expected monotonic pull increase, got ${pulls.join(', ')}`);

  // Verify no radial folding (dr'/dr > 0) even for 10x default mass
  const heavy = { x: 600, y: 400, mass: 0.02, createdAt: 0, strength: 1 };
  let prevR = -1;
  for (let r = 0; r <= 200; r += 2) {
    const pt = physics.deformSpacetimePoint(600 + r, 400, 30000, [heavy]);
    const mappedR = pt.x - 600;
    assert.ok(mappedR > prevR, `Grid folding detected at r=${r}: mappedR=${mappedR} <= prevR=${prevR}`);
    prevR = mappedR;
  }
});

test('D. Renderers: WebGL2 GLSL and Canvas 2D share causal Schwarzschild model and pass 8/8 geodesic checks', () => {
  assert.equal(physics.geodesicValidation.passed, 8, 'All 8 Schwarzschild geodesic checks must pass');
  assert.equal(physics.geodesicValidation.total, 8);

  assert.match(siteScript, /vec2 staticSchwarzschildDisplacement\(vec2 p\)/);
  assert.doesNotMatch(siteScript, /gravitationalWaveDisplacement|gravitationalWaveProjection/);
  assert.match(siteScript, /uniform float u_displayScale;/);
  assert.match(siteScript, /uniform float u_staticPull;/);
  assert.match(siteScript, /1\.0-1\.0\/\(psi\*psi\)/);
});

test('E. Regression: valid JS syntax, bumped SW cache, intact benchmark & Feedback Lab', () => {
  fs.writeFileSync('/tmp/test_site.js', siteScript);
  fs.writeFileSync('/tmp/test_feedback.js', feedbackScript);
  execFileSync(process.execPath, ['--check', '/tmp/test_site.js']);
  execFileSync(process.execPath, ['--check', '/tmp/test_feedback.js']);

  const swMatch = swJs.match(/spacetime-shell-v(\d+)/);
  assert.ok(swMatch && Number(swMatch[1]) >= 51, `Expected sw.js cache >= v51, got ${swMatch?.[0]}`);

  assert.match(siteScript, /spacetime_benchmarks/, 'Supabase benchmark endpoint preserved');
  assert.match(siteScript, /gesturestart/, 'Mac trackpad gesture handling preserved');
  assert.match(siteScript, /touchstart/, 'iPhone touch handling preserved');
  assert.match(siteScript, /resetView/, 'Reset view preserved');
  assert.match(siteScript, /visibilitychange/, 'Battery-saving visibility pause preserved');
  assert.match(siteScript, /serviceWorker\.register/, 'PWA service worker registration preserved');
});

test('F. Infinite grid: world/screen round-trip, phase snapping, and LOD scaling cover viewport at any pan/zoom/rotation', () => {
  physics.resetView();
  physics.setCam({ viewX: 4321, viewY: -2789, viewScale: 0.45, viewRotation: 37, camComX: -1200, camComY: 950 });

  for (const sx of [0, 300, 600, 1200]) {
    for (const sy of [0, 200, 400, 800]) {
      const w = physics.worldPoint(sx, sy);
      const s = physics.screenPoint(w.x, w.y);
      assert.ok(Math.abs(s.x - sx) < 1e-9, `Round-trip x mismatch: ${s.x} vs ${sx}`);
      assert.ok(Math.abs(s.y - sy) < 1e-9, `Round-trip y mismatch: ${s.y} vs ${sy}`);
    }
  }

  const cam = physics.gridCameraState();
  assert.equal(cam.lodScale, 2, 'At viewScale=0.45, LOD scale should double to 2');
  assert.ok(cam.lodScale * physics.cam.viewScale >= 0.72, 'Screen-space LOD spacing must stay >= 0.72x base');
  const center = physics.worldPoint(physics.cam.width * 0.5, physics.cam.height * 0.5);
  assert.ok(Math.abs(center.x - cam.originX) <= cam.lodStep + 1e-9, 'Grid originX snapped within 1 lodStep of viewport center');
  assert.ok(Math.abs(center.y - cam.originY) <= cam.lodStep + 1e-9, 'Grid originY snapped within 1 lodStep of viewport center');
});

test('G. N-body Newtonian gravity: binary orbit, curvature tracking, and momentum-conserving collision merge', () => {
  physics.resetView();
  const m = physics.newMassValue;
  const r = 120;
  // Circular binary orbit: v = sqrt(mu / (4 * r)) around COM
  const mu = physics.N_BODY_MU_PER_MASS * m;
  const vOrb = Math.sqrt(mu / (2 * (2 * r)));
  physics.createMass(-r, 0, m, 0, vOrb);
  physics.createMass(r, 0, m, 0, -vOrb);
  assert.equal(physics.persistentMasses.length, 2);

  // Advance 20 frames of 16ms (~0.32s)
  for (let i = 0; i < 20; i++) physics.advanceMasses(16);
  const b0 = physics.persistentMasses[0];
  const b1 = physics.persistentMasses[1];
  const sep = Math.hypot(b1.x - b0.x, b1.y - b0.y);
  assert.ok(Math.abs(sep - 2 * r) < 25, `Expected near-circular orbit separation ~${2 * r}, got ${sep.toFixed(1)}`);
  assert.ok(Math.abs(b0.y) > 5 && Math.abs(b1.y) > 5, 'Planets must move along orbit');
  assert.ok(b0.trail.length >= 2 && b1.trail.length >= 2, 'Orbital trails must be recorded');

  // Verify curvature field updates with moving planet positions
  const d0 = physics.staticSchwarzschildDisplacement(b0.x - 55, b0.y, 10000, physics.persistentMasses);
  assert.ok(Math.hypot(d0.x, d0.y) > 10, 'Curvature must track moving planet position');

  // Test momentum-conserving collision merge
  physics.resetView();
  physics.createMass(0, 0, m, 100, -40);
  physics.createMass(3, 0, 2 * m, -20, 50);
  const pxBefore = m * 100 + 2 * m * (-20);
  const pyBefore = m * (-40) + 2 * m * 50;
  physics.advanceMasses(16);
  assert.equal(physics.persistentMasses.length, 1, 'Colliding planets should merge into 1');
  const merged = physics.persistentMasses[0];
  assert.ok(Math.abs(merged.mass - 3 * m) < 1e-12, 'Mass conserved on merge');
  assert.ok(Math.abs(merged.mass * merged.vx - pxBefore) < 1e-9, 'Linear momentum X conserved on merge');
  assert.ok(Math.abs(merged.mass * merged.vy - pyBefore) < 1e-9, 'Linear momentum Y conserved on merge');
});

test('H. Slingshot launch and COM camera preservation on first planet', () => {
  physics.resetView();
  physics.setCam({ viewX: 0, viewY: 0, viewScale: 1, viewRotation: 0, camComX: 600, camComY: 400 });
  const before = physics.screenPoint(320, 240);
  physics.createMass(320, 240, physics.newMassValue, 45 * physics.SLINGSHOT_SPEED_SCALE, -30 * physics.SLINGSHOT_SPEED_SCALE);
  const after = physics.screenPoint(320, 240);
  assert.ok(Math.abs(after.x - before.x) < 1e-9 && Math.abs(after.y - before.y) < 1e-9, 'First planet creation must preserve screen position');
  assert.equal(physics.cam.camComX, 320, 'camComX initialized to first planet');
  assert.equal(physics.cam.camComY, 240, 'camComY initialized to first planet');
});

test('I. Orbital tutorial: elliptical, parabolic (with snap), hyperbolic classification and arc prediction', () => {
  physics.resetView();
  const m = physics.newMassValue;
  assert.equal(physics.classifyLaunchOrbit(100, 0, 0, 50, m), null, 'No attractor in empty space');

  physics.createMass(0, 0, m, 0, 0);
  const r = 155;
  const rSoft = Math.sqrt(r * r + 196);
  const mu = physics.N_BODY_MU_PER_MASS * (2 * m);
  const vEsc = Math.sqrt(2 * mu / rSoft);

  const ell = physics.classifyLaunchOrbit(r, 0, 0, 0.72 * vEsc, m);
  assert.equal(ell.regime, 'elliptical');
  assert.ok(ell.e < 1, `Expected elliptical e < 1, got ${ell.e}`);

  // Near escape speed (within 4.5% band) snaps gently to parabolic
  const par = physics.classifyLaunchOrbit(r, 0, 0, 0.98 * vEsc, m);
  assert.equal(par.regime, 'parabolic');
  assert.equal(par.snapped, true, 'Expected parabolic magnetic snap near v_esc');
  assert.ok(Math.abs(par.ratio - 1) < 1e-9, 'Snapped ratio should be 1.0');

  const hyp = physics.classifyLaunchOrbit(r, 0, 0, 1.38 * vEsc, m);
  assert.equal(hyp.regime, 'hyperbolic');
  assert.ok(hyp.e > 1, `Expected hyperbolic e > 1, got ${hyp.e}`);

  // Radial infall below escape speed should still be classified as bound/elliptical
  const radialBound = physics.classifyLaunchOrbit(r, 0, -0.5 * vEsc, 0, m);
  assert.equal(radialBound.regime, 'elliptical');

  assert.ok(physics.predictLaunchArc(r, 0, 0, 0.72 * vEsc, m) > 20, 'Predicted launch arc should sample trajectory points');
});
