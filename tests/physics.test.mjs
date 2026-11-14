import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { execFileSync } from 'node:child_process';

const rootDir = path.resolve(import.meta.dirname, '..');
const indexHtml = fs.readFileSync(path.join(rootDir, 'index.html'), 'utf8');
const feedbackHtml = fs.readFileSync(path.join(rootDir, 'feedback-lab', 'index.html'), 'utf8');
const spinWavesHtml = fs.readFileSync(path.join(rootDir, 'spin-waves', 'index.html'), 'utf8');
const swJs = fs.readFileSync(path.join(rootDir, 'sw.js'), 'utf8');

function extractInlineScript(html, label) {
  const matches = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)];
  assert.equal(matches.length, 1, `Expected 1 inline script in ${label}, found ${matches.length}`);
  return matches[0][1];
}

const siteScript = extractInlineScript(indexHtml, 'index.html');
const feedbackScript = extractInlineScript(feedbackHtml, 'feedback-lab/index.html');
const spinWavesScript = extractInlineScript(spinWavesHtml, 'spin-waves/index.html');

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
  fs.writeFileSync('/tmp/test_spinwaves.js', spinWavesScript);
  execFileSync(process.execPath, ['--check', '/tmp/test_site.js']);
  execFileSync(process.execPath, ['--check', '/tmp/test_feedback.js']);
  execFileSync(process.execPath, ['--check', '/tmp/test_spinwaves.js']);

  const swMatch = swJs.match(/spacetime-shell-v(\d+)/);
  assert.ok(swMatch && Number(swMatch[1]) >= 53, `Expected sw.js cache >= v53, got ${swMatch?.[0]}`);
  assert.match(indexHtml, /\.\/spin-waves\//, 'Spin Waves navigation button preserved');

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

function loadSpinWavesContext() {
  const noop = () => {};
  const makeCtx = () => new Proxy({
    createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) })
  }, {
    get(target, prop) { return prop in target ? target[prop] : noop; }
  });
  const makeEl = () => ({
    style: {},
    classList: { add: noop, remove: noop, toggle: noop },
    setAttribute: noop,
    addEventListener: noop,
    getContext: () => makeCtx(),
    querySelector: () => makeEl(),
    textContent: '',
    innerHTML: '',
    value: 0,
    width: 280,
    height: 54
  });
  const sandbox = {
    console,
    Math,
    Number,
    Array,
    Float32Array,
    Uint8Array,
    Uint8ClampedArray,
    innerWidth: 1200,
    innerHeight: 800,
    devicePixelRatio: 2,
    localStorage: { getItem: () => null, setItem: noop },
    requestAnimationFrame: () => 0,
    addEventListener: noop,
    document: {
      getElementById: () => makeEl(),
      createElement: () => makeEl(),
      querySelectorAll: () => [],
      addEventListener: noop
    },
    window: null
  };
  sandbox.window = sandbox;

  const instrumented = spinWavesScript.replace(
    /\}\)\(\);\s*$/,
    `  globalThis.__spin = {
    NX, NY, N, J_EX, DEMAG_NZ, FREQ_TO_RAD,
    mx, my, mz, activeMask, modeTemplates, MODE_SPECS,
    buildMaskAndModes, seedSpins, exciteEigenmode, stepLlg,
    modeFreqGHz, parametricThreshold, besselJ,
    setState(patch) {
      if ('geometryMode' in patch) geometryMode = patch.geometryMode;
      if ('driveMode' in patch) driveMode = patch.driveMode;
      if ('driveFreq' in patch) driveFreq = patch.driveFreq;
      if ('driveAmp' in patch) driveAmp = patch.driveAmp;
      if ('biasH0' in patch) biasH0 = patch.biasH0;
      if ('damping' in patch) damping = patch.damping;
      if ('simTime' in patch) simTime = patch.simTime;
    },
    getState() { return { geometryMode, driveMode, driveFreq, driveAmp, biasH0, damping, simTime }; },
    meanTransverseAmp() {
      let sum = 0, cnt = 0;
      for (let i = 0; i < N; i++) {
        if (!activeMask[i]) continue;
        sum += Math.hypot(my[i], mz[i]); cnt++;
      }
      return cnt ? sum / cnt : 0;
    },
    modeProjection(mIdx) {
      let dy = 0, dz = 0;
      const tpl = modeTemplates[mIdx];
      for (let i = 0; i < N; i++) {
        if (!activeMask[i]) continue;
        dy += my[i] * tpl[i]; dz += mz[i] * tpl[i];
      }
      return Math.hypot(dy, dz);
    }
  };
})();`
  );
  const context = vm.createContext(sandbox);
  vm.runInContext(instrumented, context);
  return context.__spin;
}

const spin = loadSpinWavesContext();

test('J. Spin-wave model: unit spin norm conservation and eigenmode orthonormality', () => {
  for (const geom of ['disk', 'rect']) {
    spin.setState({ geometryMode: geom });
    spin.buildMaskAndModes();
    for (let a = 0; a < spin.MODE_SPECS.length; a++) {
      let normSq = 0;
      for (let i = 0; i < spin.N; i++) normSq += spin.modeTemplates[a][i] ** 2;
      assert.ok(Math.abs(normSq - 1) < 1e-5, `${geom} mode ${a} not normalized: ${normSq}`);
      for (let b = a + 1; b < spin.MODE_SPECS.length; b++) {
        let dot = 0;
        for (let i = 0; i < spin.N; i++) dot += spin.modeTemplates[a][i] * spin.modeTemplates[b][i];
        assert.ok(Math.abs(dot) < 0.06, `${geom} modes ${a},${b} not orthogonal: dot=${dot.toFixed(4)}`);
      }
    }
  }

  spin.setState({ geometryMode: 'disk', driveMode: 'parallel', driveAmp: 0.35, damping: 0.04 });
  spin.buildMaskAndModes();
  spin.seedSpins(0.03);
  for (let s = 0; s < 300; s++) spin.stepLlg(0.01);
  for (let i = 0; i < spin.N; i++) {
    if (!spin.activeMask[i]) continue;
    const len = Math.hypot(spin.mx[i], spin.my[i], spin.mz[i]);
    assert.ok(Math.abs(len - 1) < 1e-5, `Spin ${i} norm drifted to ${len}`);
  }
});

test('K. Spin-wave model: eigenmode dispersion & mode purity (simulated vs analytical frequency)', () => {
  spin.setState({ geometryMode: 'disk', biasH0: 0.45, damping: 0.0, driveAmp: 0.0, simTime: 0 });
  spin.buildMaskAndModes();

  for (const mIdx of [0, 1, 3]) {
    spin.exciteEigenmode(mIdx);
    spin.setState({ damping: 0.0, driveAmp: 0.0, simTime: 0 });
    const fExpected = spin.modeFreqGHz(mIdx);

    const proj = () => {
      let sum = 0;
      for (let i = 0; i < spin.N; i++) if (spin.activeMask[i]) sum += spin.my[i] * spin.modeTemplates[mIdx][i];
      return sum;
    };
    const crossings = [];
    let prev = proj();
    const dt = 0.01;
    for (let step = 1; step <= 3000; step++) {
      spin.stepLlg(dt);
      const cur = proj();
      if (prev <= 0 && cur > 0) {
        const frac = -prev / (cur - prev);
        crossings.push((step - 1 + frac) * dt);
        if (crossings.length >= 4) break;
      }
      prev = cur;
    }
    assert.ok(crossings.length >= 3, `Mode ${mIdx} did not complete enough oscillations`);
    const avgPeriod = (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1);
    const omegaSim = (2 * Math.PI) / avgPeriod;
    const fSim = omegaSim / spin.FREQ_TO_RAD;
    const relErr = Math.abs(fSim - fExpected) / fExpected;
    assert.ok(relErr < 0.065, `Mode ${mIdx} frequency mismatch: sim=${fSim.toFixed(3)} vs expected=${fExpected.toFixed(3)} GHz`);

    const pSelf = spin.modeProjection(mIdx);
    let totalSq = 0;
    for (let k = 0; k < spin.MODE_SPECS.length; k++) totalSq += spin.modeProjection(k) ** 2;
    assert.ok(pSelf / Math.sqrt(totalSq) > 0.95, `Mode ${mIdx} leaked into other modes: purity=${(pSelf / Math.sqrt(totalSq)).toFixed(3)}`);
  }
});

test('L. Spin-wave model: Gilbert damping energy dissipation when drive is removed', () => {
  spin.setState({ geometryMode: 'disk', biasH0: 0.45, driveAmp: 0.0 });
  spin.buildMaskAndModes();

  spin.exciteEigenmode(1);
  spin.setState({ damping: 0.0, driveAmp: 0.0 });
  const a0Zero = spin.meanTransverseAmp();
  for (let s = 0; s < 1200; s++) spin.stepLlg(0.01);
  const a1Zero = spin.meanTransverseAmp();
  assert.ok(a1Zero >= a0Zero * 0.42, `Zero-damping amplitude collapsed unexpectedly: ${a1Zero.toFixed(4)} vs ${a0Zero.toFixed(4)}`);

  spin.exciteEigenmode(1);
  spin.setState({ damping: 0.04, driveAmp: 0.0 });
  for (let s = 0; s < 2400; s++) spin.stepLlg(0.01);
  const a1Damp = spin.meanTransverseAmp();
  assert.ok(a1Damp < a0Zero * 0.35, `Expected strong zero-drive decay at alpha=0.04: ${a1Damp.toFixed(4)} vs ${a0Zero.toFixed(4)}`);
});

test('M. Spin-wave model: parametric excitation threshold (decay below h_th, exponential growth above h_th)', () => {
  spin.setState({ geometryMode: 'disk', driveMode: 'parallel', biasH0: 0.45, damping: 0.040 });
  spin.buildMaskAndModes();
  const hTh = spin.parametricThreshold();
  const fPump = 2 * spin.modeFreqGHz(3);

  // Below threshold: h_p = 0.25 * h_th -> resonant mode 3 must decay
  spin.setState({ driveFreq: fPump, driveAmp: 0.25 * hTh, simTime: 0 });
  spin.seedSpins(0.03);
  const subStart = spin.modeProjection(3);
  for (let s = 0; s < 1800; s++) spin.stepLlg(0.01);
  const subEnd = spin.modeProjection(3);
  assert.ok(subEnd < subStart * 0.55, `Below threshold (h_p=${(0.25 * hTh).toFixed(3)} < h_th=${hTh.toFixed(3)}), mode 3 should decay: ${subStart.toFixed(4)} -> ${subEnd.toFixed(4)}`);

  // Above threshold: h_p = 2.5 * h_th -> parametric instability must amplify resonant mode 3!
  spin.setState({ driveFreq: fPump, driveAmp: 2.5 * hTh, simTime: 0 });
  spin.seedSpins(0.03);
  const supStart = spin.modeProjection(3);
  for (let s = 0; s < 3000; s++) spin.stepLlg(0.01);
  const supEnd = spin.modeProjection(3);
  assert.ok(supEnd > supStart * 2.2, `Above threshold (h_p=${(2.5 * hTh).toFixed(3)} > h_th=${hTh.toFixed(3)}), mode 3 should grow: ${supStart.toFixed(4)} -> ${supEnd.toFixed(4)}`);
});
