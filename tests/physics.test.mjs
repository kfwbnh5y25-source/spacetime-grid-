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

// Load the pure physics & deformation functions from index.html in an isolated VM
// with minimal DOM/browser stubs so tests execute the actual production code.
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
      addEventListener: noop,
      body: { appendChild: noop }
    },
    window: null
  };
  sandbox.window = sandbox;

  // Expose internal physics functions at the end of the IIFE for inspection
  const instrumented = siteScript.replace(
    /\}\)\(\);\s*$/,
    `  globalThis.__physics = {
    DISPLAY_SCHWARZSCHILD_SCALE,
    STATIC_WELL_PULL_PX,
    STATIC_WELL_SOFTENING_PX,
    staticSchwarzschildDisplacement,
    gravitationalWaveProjection,
    deformSpacetimePoint,
    runGeodesicValidation,
    geodesicValidation,
    newMassValue,
    gravityG
  };
})();`
  );

  const context = vm.createContext(sandbox);
  vm.runInContext(instrumented, context);
  return context.__physics;
}

const physics = loadPhysicsContext();

test('A. Static curvature: visible, symmetric, smooth, and time-independent', () => {
  const mass = { x: 600, y: 400, mass: physics.newMassValue, createdAt: 0, strength: 1 };
  const timestamps = [0, 500, 2000, 10000, 60000];
  const samplePoints = [
    { x: 640, y: 400 },
    { x: 560, y: 400 },
    { x: 600, y: 440 },
    { x: 600, y: 360 },
    { x: 680, y: 400 },
    { x: 760, y: 400 }
  ];

  // Static displacement must be identical across all animation timestamps
  const baseline = samplePoints.map((pt) => physics.staticSchwarzschildDisplacement(pt.x, pt.y, [mass]));
  for (const t of timestamps) {
    samplePoints.forEach((pt, idx) => {
      const d = physics.staticSchwarzschildDisplacement(pt.x, pt.y, [mass], t);
      assert.equal(d.x, baseline[idx].x, `Static dx changed at t=${t}`);
      assert.equal(d.y, baseline[idx].y, `Static dy changed at t=${t}`);
    });
  }

  // Radial symmetry around the stationary mass (40 px in +x, -x, +y, -y)
  const dPx = physics.staticSchwarzschildDisplacement(640, 400, [mass]);
  const dNx = physics.staticSchwarzschildDisplacement(560, 400, [mass]);
  const dPy = physics.staticSchwarzschildDisplacement(600, 440, [mass]);
  const dNy = physics.staticSchwarzschildDisplacement(600, 360, [mass]);
  assert.ok(dPx.x < 0 && Math.abs(dPx.y) < 1e-12, 'Pull on +x axis points toward mass (-x)');
  assert.ok(Math.abs(dPx.x + dNx.x) < 1e-12, 'Symmetric between +x and -x');
  assert.ok(Math.abs(dPx.x - dPy.y) < 1e-12, 'Symmetric between x and y axes');
  assert.ok(Math.abs(dPy.y + dNy.y) < 1e-12, 'Symmetric between +y and -y');

  // Unmistakable visibility (~12-18 px peak pull near mass, decaying smoothly with distance)
  const pull40 = Math.hypot(dPx.x, dPx.y);
  const pull120 = Math.hypot(...Object.values(physics.staticSchwarzschildDisplacement(720, 400, [mass])));
  const pull300 = Math.hypot(...Object.values(physics.staticSchwarzschildDisplacement(900, 400, [mass])));
  assert.ok(pull40 >= 10 && pull40 <= 20, `Expected visible 10-20 px pull at 40 px, got ${pull40.toFixed(2)} px`);
  assert.ok(pull120 > 5 && pull120 < pull40, `Expected smooth decay at 120 px, got ${pull120.toFixed(2)} px`);
  assert.ok(pull300 < pull120 && pull300 > 1, `Expected continued 1/r decay at 300 px, got ${pull300.toFixed(2)} px`);
});

test('B. Gravitational wave: causal propagation, 2D + polarization, finite-time decay, static curvature persists', () => {
  const mass = { x: 600, y: 400, mass: physics.newMassValue, createdAt: 0, strength: 1 };

  // Ahead of the light cone (r = 360 px, c = 180 px/s => arrival at t = 2.0 s)
  const preArrival = physics.gravitationalWaveProjection(960, 400, 1500, [mass]);
  assert.equal(preArrival.x, 0, 'Wave must be zero ahead of light cone');
  assert.equal(preArrival.y, 0, 'Wave must be zero ahead of light cone');

  // During packet passage at r = 180 px (arrival at 1.0 s, packet peaks around retarded ~ 0.85 s => t ~ 1.85-2.1 s)
  let peakX = 0;
  let peakY = 0;
  for (let tMs = 1200; tMs <= 2600; tMs += 50) {
    const wx = physics.gravitationalWaveProjection(780, 400, tMs, [mass]);
    const wy = physics.gravitationalWaveProjection(600, 580, tMs, [mass]);
    peakX = Math.max(peakX, Math.abs(wx.x));
    peakY = Math.max(peakY, Math.abs(wy.y));
    // + polarization symmetry: at equal t, +x and +y have opposite signs
    assert.ok(wx.x * wy.y <= 1e-12, `Expected opposite signs for + polarization at t=${tMs}`);
  }
  assert.ok(peakX > 2 && peakY > 2, `Expected visible 2D wave amplitude, got peakX=${peakX.toFixed(2)}, peakY=${peakY.toFixed(2)}`);

  // After packet passes (t = 15000 ms => retarded = 14 s at r = 180 px)
  const lateWave = physics.gravitationalWaveProjection(780, 400, 15000, [mass]);
  assert.ok(Math.hypot(lateWave.x, lateWave.y) < 1e-6, 'Wave must decay to zero after passing');

  // Total deformation after wave disappears equals static Schwarzschild curvature
  const totalLate = physics.deformSpacetimePoint(640, 400, 15000, [mass]);
  const staticOnly = physics.staticSchwarzschildDisplacement(640, 400, [mass]);
  assert.ok(Math.abs((totalLate.x - 640) - staticOnly.x) < 1e-6, 'Static curvature must remain after wave passes');
  assert.ok(Math.abs((totalLate.y - 400) - staticOnly.y) < 1e-6, 'Static curvature must remain after wave passes');
});

test('C. Mass scaling: monotonic curvature increase across 3 mass strengths and no line folding', () => {
  const masses = [0.001, 0.002, 0.006];
  const pulls = masses.map((m) => {
    const src = { x: 600, y: 400, mass: m, createdAt: 0, strength: 1 };
    const d = physics.staticSchwarzschildDisplacement(640, 400, [src]);
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

test('D. Renderers: WebGL2 GLSL and Canvas 2D share identical model and pass 8/8 geodesic checks', () => {
  assert.equal(physics.geodesicValidation.passed, 8, 'All 8 Schwarzschild geodesic checks must pass');
  assert.equal(physics.geodesicValidation.total, 8);

  // Verify WebGL2 shader contains the matching static Schwarzschild & 2D + polarization wave functions
  assert.match(siteScript, /vec2 staticSchwarzschildDisplacement\(vec2 p\)/);
  assert.match(siteScript, /vec2 gravitationalWaveDisplacement\(vec2 p\)/);
  assert.match(siteScript, /uniform float u_displayScale;/);
  assert.match(siteScript, /uniform float u_staticPull;/);
  assert.match(siteScript, /1\.0-1\.0\/\(psi\*psi\)/);
  assert.match(siteScript, /vec2\(d\.x,-d\.y\)\*reg/);
});

test('E. Regression: valid JS syntax, bumped SW cache, intact benchmark & Feedback Lab', () => {
  // Syntax check via node --check on extracted scripts
  fs.writeFileSync('/tmp/test_site.js', siteScript);
  fs.writeFileSync('/tmp/test_feedback.js', feedbackScript);
  execFileSync(process.execPath, ['--check', '/tmp/test_site.js']);
  execFileSync(process.execPath, ['--check', '/tmp/test_feedback.js']);

  // Service worker cache bumped past v46
  const swMatch = swJs.match(/spacetime-shell-v(\d+)/);
  assert.ok(swMatch && Number(swMatch[1]) >= 47, `Expected sw.js cache >= v47, got ${swMatch?.[0]}`);

  // Regression hooks preserved
  assert.match(siteScript, /spacetime_benchmarks/, 'Supabase benchmark endpoint preserved');
  assert.match(siteScript, /gesturestart/, 'Mac trackpad gesture handling preserved');
  assert.match(siteScript, /touchstart/, 'iPhone touch handling preserved');
  assert.match(siteScript, /resetView/, 'Reset view preserved');
  assert.match(siteScript, /visibilitychange/, 'Battery-saving visibility pause preserved');
  assert.match(siteScript, /serviceWorker\.register/, 'PWA service worker registration preserved');
});
