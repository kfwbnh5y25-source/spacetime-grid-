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
    querySelectorAll: () => [],
    textContent: '',
    innerHTML: '',
    value: 0,
    disabled: false,
    width: 280,
    height: 54
  });
  const sandbox = {
    console,
    Math,
    Number,
    Array,
    Float32Array,
    Float64Array,
    Int8Array,
    Int16Array,
    Int32Array,
    Uint8Array,
    Uint16Array,
    Uint8ClampedArray,
    innerWidth: 1200,
    innerHeight: 800,
    devicePixelRatio: 2,
    localStorage: { getItem: () => null, setItem: noop },
    requestAnimationFrame: () => 0,
    addEventListener: noop,
    atob: s => Buffer.from(s, 'base64').toString('binary'),
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
    M, NX, NY, N, BDC_MT, MU0_MS_MT, F_M, FREQ_TO_RAD, SPEC_LEN, SPEC_DT,
    TUKEY_RAMP_NS, TUKEY_RAMP_TAU, THETA_PER_K,
    NACT, actPix, mx, my, mz, m0x, m0y, m0z, m0xFull, m0yFull, m0zFull, activeMask,
    phiReX, phiReY, phiReZ, phiImX, phiImY, phiImZ,
    aRe, aIm, wBarArr, Kp, Km, N4A, P4B, terms, isRwaFast, setNmmMode,
    MODE_SPECS, spectrum, reconstructSpatialSpins,
    thermalTheta, thermalRmsAmp, reseedRng,
    nsToTau, usToTau, simTimeNs, simTimeUs,
    tukeyRamp, tone1EnvAt, tone2EnvAt,
    seedSpins, exciteEigenmode, pluckAt, stepLlg, updateSpectrum,
    modeFreqGHz, modeGammaGHz, modeGamma, parametricThreshold, parametricThresholdMT, arnoldThreshold,
    cycleAvgModeAmp, modeRmsAmp,
    selfNfsCoeff, mutualNfsCoeff, mutualNfsCoeffGHz, modeNfsShiftGHz, modeMutualNfsShiftGHz,
    steadyStateAmp, analyzeTwoMode, steadyState2Mode, fastForwardSteadyState, fastForwardStepsPerFrame,
    applyPreset, swapTonesAction, SPEED_FRAME_STEPS,
    swLog, recordLogEvent, recordTelemetrySample, buildEvaluationReport,
    drawModeTimePlot, updatePlotWindow, recordPlotSample,
    getPlotStats: () => ({ count: plotHistCount, winTau: plotWinTau, originTau: plotOriginTau, lastPlotTau }),
    setState(patch) {
      if ('driveMode' in patch) driveMode = patch.driveMode;
      if ('driveFreq' in patch) driveFreq = patch.driveFreq;
      if ('driveAmp' in patch) driveAmp = patch.driveAmp;
      if ('damping' in patch) damping = patch.damping;
      if ('dualTone' in patch) dualTone = patch.dualTone;
      if ('driveFreq2' in patch) driveFreq2 = patch.driveFreq2;
      if ('driveAmp2' in patch) driveAmp2 = patch.driveAmp2;
      if ('tone2DelayNs' in patch) tone2DelayNs = patch.tone2DelayNs;
      if ('tone2WaitSS' in patch) {
        tone2WaitSS = patch.tone2WaitSS;
        tone2Released = false; tone2ReleaseT = Infinity; t1WaitStable = 0; lastWaitCheckT = simTime;
      }
      if ('drivePhase1' in patch) drivePhase1 = patch.drivePhase1;
      if ('drivePhase2' in patch) drivePhase2 = patch.drivePhase2;
      if ('temperatureK' in patch) temperatureK = patch.temperatureK;
      if ('thermalSeed' in patch) reseedRng(patch.thermalSeed);
      if ('arnoldTarget' in patch) arnoldTarget = patch.arnoldTarget;
      if ('playbackSpeed' in patch) playbackSpeed = patch.playbackSpeed;
      if ('simTime' in patch) { simTime = patch.simTime; resetModeHistory(); recordModeSample(); }
    },
    getState() {
      return {
        driveMode, driveFreq, driveAmp, damping, simTime, dualTone, driveFreq2, driveAmp2,
        tone2DelayNs, tone2WaitSS, tone2Released, tone2ReleaseT,
        drivePhase1, drivePhase2, temperatureK, playbackSpeed, isRwa: isRwaFast()
      };
    },
    meanTransverseAmp() {
      let sum = 0, cnt = 0;
      for (let i = 0; i < N; i++) {
        if (!activeMask[i]) continue;
        sum += Math.hypot(my[i], mz[i]); cnt++;
      }
      return cnt ? sum / cnt : 0;
    },
    symplecticInnerProduct(mIdx) {
      let re = 0, im = 0;
      const off = mIdx * NACT;
      const invN = 1 / NACT;
      for (let p = 0; p < NACT; p++) {
        const pix = actPix[p];
        const dmx = mx[pix] - m0x[p], dmy = my[pix] - m0y[p], dmz = mz[pix] - m0z[p];
        const cx = m0y[p] * dmz - m0z[p] * dmy;
        const cy = m0z[p] * dmx - m0x[p] * dmz;
        const cz = m0x[p] * dmy - m0y[p] * dmx;
        re += -(phiImX[off + p] * cx + phiImY[off + p] * cy + phiImZ[off + p] * cz);
        im += -(phiReX[off + p] * cx + phiReY[off + p] * cy + phiReZ[off + p] * cz);
      }
      return Math.hypot(wBarArr[mIdx] * re * invN, wBarArr[mIdx] * im * invN);
    }
  };
})();`
  );
  const context = vm.createContext(sandbox);
  vm.runInContext(instrumented, context);
  return context.__spin;
}

const spin = loadSpinWavesContext();

test("J. Spin-wave NMM: 29 modes (>1.6 GHz), uniform field only, unit spin norm & symplectic B0 reconstruction", () => {
  assert.equal(spin.M, 29, 'Must include strictly 29 normal modes above 1.6 GHz');
  assert.equal(spin.MODE_SPECS.length, 29);
  for (let m = 0; m < spin.M; m++) {
    assert.ok(spin.MODE_SPECS[m].fGHz > 1.6, `Mode ${m} frequency ${spin.MODE_SPECS[m].fGHz} GHz must exceed 1.6 GHz`);
  }
  assert.ok(spin.MODE_SPECS[0].fGHz > 1.70 && spin.MODE_SPECS[28].fGHz < 2.08);
  assert.doesNotMatch(spinWavesHtml, /id="slChi|id="nikBtnZeroMut"|buildMaskAndModes|twoModeIsolate/,
    'Legacy custom sliders/gradients/shortcuts must be completely removed');

  // Symplectic B0 round-trip on Kittel (m27) and bulk (m13) modes + unit norm conservation
  spin.setState({ temperatureK: 0, driveAmp: 0, dualTone: false, driveAmp2: 0, damping: 0.002 });
  for (const mIdx of [13, 27]) {
    spin.exciteEigenmode(mIdx, 0.04);
    const proj = spin.symplecticInnerProduct(mIdx);
    assert.ok(Math.abs(proj - 0.04) / 0.04 < 0.05, `Symplectic roundtrip on m${mIdx}: expected 0.04, got ${proj.toFixed(5)}`);
  }

  for (const full of [false, true]) {
    spin.setNmmMode(full);
    spin.setState({ temperatureK: 0, driveMode: 'parallel', driveFreq: 2 * spin.modeFreqGHz(13), driveAmp: 1.8 * spin.parametricThreshold(13), damping: 0.002, simTime: 0 });
    spin.seedSpins(0.01);
    for (let s = 0; s < 300; s++) spin.stepLlg(0.65, true);
    spin.reconstructSpatialSpins();
    for (let p = 0; p < spin.NACT; p++) {
      const pix = spin.actPix[p];
      const len = Math.hypot(spin.mx[pix], spin.my[pix], spin.mz[pix]);
      assert.ok(Math.abs(len - 1) < 1e-5, `Spin pixel ${pix} norm drifted to ${len} (full=${full})`);
    }
  }
  spin.setNmmMode(false);
});

test('K. Spin-wave NMM: exact linear eigenmode frequencies and modal purity', () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, damping: 0.0, driveAmp: 0.0, dualTone: false, driveAmp2: 0.0, simTime: 0 });
  for (const mIdx of [10, 13, 27]) {
    spin.exciteEigenmode(mIdx, 0.005);
    spin.setState({ damping: 0.0, driveAmp: 0.0, simTime: 0 });
    const fExpected = spin.modeFreqGHz(mIdx);
    const crossings = [];
    let prev = spin.aRe[mIdx];
    const dt = 0.25;
    for (let step = 1; step <= 800; step++) {
      spin.stepLlg(dt, true);
      const cur = spin.aRe[mIdx];
      if (prev <= 0 && cur > 0) {
        const frac = -prev / (cur - prev);
        crossings.push((step - 1 + frac) * dt);
        if (crossings.length >= 6) break;
      }
      prev = cur;
    }
    assert.ok(crossings.length >= 4, `Mode ${mIdx} did not complete enough linear cycles`);
    const avgPeriod = (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1);
    const fSim = ((2 * Math.PI) / avgPeriod) / spin.FREQ_TO_RAD;
    const relErr = Math.abs(fSim - fExpected) / fExpected;
    assert.ok(relErr < 0.02, `m${mIdx} freq mismatch: sim=${fSim.toFixed(4)} vs ${fExpected.toFixed(4)} GHz`);
    let otherSq = 0;
    for (let k = 0; k < spin.M; k++) if (k !== mIdx) otherSq += spin.aRe[k] ** 2 + spin.aIm[k] ** 2;
    assert.ok(Math.sqrt(otherSq) < 1e-6, `m${mIdx} leaked into other modes in linear zero-drive precession`);
  }
});

test('L. Spin-wave NMM: Gilbert damping energy dissipation when drive is removed', () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, dualTone: false, driveAmp2: 0 });
  spin.exciteEigenmode(13, 0.04);
  spin.setState({ damping: 0.0, driveAmp: 0.0, simTime: 0 });
  const a0 = spin.cycleAvgModeAmp(13);
  for (let s = 0; s < 1500; s++) spin.stepLlg(0.65, true);
  assert.ok(spin.cycleAvgModeAmp(13) > a0 * 0.99, 'Zero-damping amplitude should be conserved in RWA core');

  spin.exciteEigenmode(13, 0.04);
  spin.setState({ damping: 0.004, driveAmp: 0.0, simTime: 0 });
  for (let s = 0; s < 3000; s++) spin.stepLlg(0.65, true); // 1950 tau
  const expectedDecay = Math.exp(-spin.modeGamma(13) * 1950);
  const actualDecay = spin.cycleAvgModeAmp(13) / a0;
  assert.ok(actualDecay < 0.35 && Math.abs(actualDecay - expectedDecay) < 0.04,
    `Gilbert decay mismatch: actual=${actualDecay.toFixed(4)} vs exp(-Gamma*tau)=${expectedDecay.toFixed(4)}`);
});

test('M. Spin-wave NMM: parallel parametric excitation threshold (decay below h_th, growth above h_th)', () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, driveMode: 'parallel', damping: 0.002, dualTone: false, driveAmp2: 0 });
  const hTh = spin.parametricThreshold(13);
  const fPump = 2 * spin.modeFreqGHz(13);

  spin.setState({ driveFreq: fPump, driveAmp: 0.4 * hTh, simTime: 0 });
  spin.seedSpins(0.005, [13]);
  const aStartSub = spin.cycleAvgModeAmp(13);
  for (let s = 0; s < 4000; s++) spin.stepLlg(0.65, true);
  assert.ok(spin.cycleAvgModeAmp(13) < aStartSub * 0.3, 'Below h_th, m13 must decay');

  spin.setState({ driveFreq: fPump, driveAmp: 2.0 * hTh, simTime: 0 });
  spin.seedSpins(0.002, [13]);
  for (let s = 0; s < 8000; s++) spin.stepLlg(0.65, true);
  assert.ok(spin.cycleAvgModeAmp(13) > 0.025, `Above h_th, m13 must amplify: got ${spin.cycleAvgModeAmp(13).toFixed(4)}`);
});

test("N. Spin-wave NMM: d'Aquino Eq. (33) steady-state amplitudes (hardening m13, softening m27) & Arnold boundary", () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, driveMode: 'parallel', damping: 0.002, dualTone: false, driveAmp2: 0, tone2DelayNs: 0, tone2WaitSS: false });

  // Hardening bulk mode m13 (N_13,13 > 0)
  const n13 = spin.mutualNfsCoeffGHz(13, 13);
  const n27 = spin.mutualNfsCoeffGHz(27, 27);
  assert.ok(n13 > 2.5, `Expected hardening N_13,13 > +2.5 GHz, got ${n13}`);
  assert.ok(n27 < -2.0, `Expected softening N_27,27 < -2.0 GHz, got ${n27}`);

  const fp13 = 2 * spin.modeFreqGHz(13) + 0.002;
  const hp13 = 1.6 * spin.parametricThreshold(13);
  spin.setState({ driveFreq: fp13, driveAmp: hp13, simTime: 0 });
  spin.seedSpins(0.005, [13]);
  spin.fastForwardSteadyState();
  const th13 = spin.steadyStateAmp(13, fp13, hp13, 0, true);
  const sim13 = spin.cycleAvgModeAmp(13);
  assert.ok(th13 > 0.04 && Math.abs(sim13 - th13) / th13 < 0.05,
    `Hardening m13 Eq.(33) mismatch: sim=${sim13.toFixed(4)} vs theory=${th13.toFixed(4)}`);

  // Softening Kittel mode m27 (N_27,27 < 0, pump slightly below 2*f27)
  const fp27 = 2 * spin.modeFreqGHz(27) - 0.002;
  const hp27 = 1.7 * spin.parametricThreshold(27);
  spin.setState({ driveFreq: fp27, driveAmp: hp27, simTime: 0 });
  spin.seedSpins(0.005, [27]);
  spin.fastForwardSteadyState();
  const th27 = spin.steadyStateAmp(27, fp27, hp27, 0, true);
  const sim27 = spin.cycleAvgModeAmp(27);
  assert.ok(th27 > 0.03 && Math.abs(sim27 - th27) / th27 < 0.05,
    `Softening m27 Eq.(33) mismatch: sim=${sim27.toFixed(4)} vs theory=${th27.toFixed(4)}`);

  // Dynamic Arnold tongue check: outside detuned tongue must decay
  const hTip = spin.parametricThreshold(13);
  assert.ok(spin.arnoldThreshold(13, 2 * spin.modeFreqGHz(13) + 0.04) > hTip * 2);
  spin.setState({ driveFreq: 2 * spin.modeFreqGHz(13) + 0.06, driveAmp: 1.5 * hTip, simTime: 0 });
  spin.seedSpins(0.005, [13]);
  for (let s = 0; s < 8000; s++) spin.stepLlg(0.65, true);
  assert.ok(spin.cycleAvgModeAmp(13) < 0.001, 'Outside Arnold tongue must decay to equilibrium');
});

test("O. Spin-wave NMM: d'Aquino Scenario 1 (S_det > 0) coupled coexistence on (13, 23) and opposite-sign (13, 27)", () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0 });

  // Weakly coupled hardening pair (13, 23)
  spin.applyPreset('coex');
  let st = spin.getState();
  const anaCoex = spin.analyzeTwoMode(13, 23, st.driveFreq, st.driveAmp, st.driveFreq2, st.driveAmp2, true);
  assert.equal(anaCoex.scenario, 1);
  assert.ok(anaCoex.S_det > 0 && anaCoex.branches['c++']?.stable);
  spin.fastForwardSteadyState();
  const s13 = spin.cycleAvgModeAmp(13), s23 = spin.cycleAvgModeAmp(23);
  assert.ok(s13 > 0.025 && s23 > 0.08, `Expected (13,23) c++ coexistence, got (${s13.toFixed(4)}, ${s23.toFixed(4)})`);
  assert.ok(Math.abs(s13 - anaCoex.pred[0]) / anaCoex.pred[0] < 0.09, `m13 mismatch: ${s13.toFixed(4)} vs ${anaCoex.pred[0].toFixed(4)}`);
  assert.ok(Math.abs(s23 - anaCoex.pred[1]) / anaCoex.pred[1] < 0.05, `m23 mismatch: ${s23.toFixed(4)} vs ${anaCoex.pred[1].toFixed(4)}`);

  // Opposite-sign pair (13 hardening, 27 softening, Fig. 5)
  spin.applyPreset('opp');
  st = spin.getState();
  const anaOpp = spin.analyzeTwoMode(13, 27, st.driveFreq, st.driveAmp, st.driveFreq2, st.driveAmp2, true);
  assert.equal(anaOpp.scenario, 1);
  spin.fastForwardSteadyState();
  const o13 = spin.cycleAvgModeAmp(13), o27 = spin.cycleAvgModeAmp(27);
  assert.ok(o13 > 0.025 && o27 > 0.035, `Expected (13,27) c++ coexistence, got (${o13.toFixed(4)}, ${o27.toFixed(4)})`);
  assert.ok(Math.abs(o13 - anaOpp.pred[0]) / anaOpp.pred[0] < 0.08);
  assert.ok(Math.abs(o27 - anaOpp.pred[1]) / anaOpp.pred[1] < 0.08);
});

test("P. Spin-wave NMM: d'Aquino Scenario 2 (S_det < 0) 3.0 µs delayed nonreciprocity, Wait T1 SS & Swap T1<->T2", () => {
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0 });
  spin.applyPreset('nonrecip');

  // Verify analytical Scenario 2 & unstable c++ branch
  const st = spin.getState();
  const ana = spin.analyzeTwoMode(13, 16, st.driveFreq, st.driveAmp, st.driveFreq2, st.driveAmp2, true);
  assert.equal(ana.scenario, 2);
  assert.ok(ana.S_det < 0, `Expected S_det < 0 on (13, 16), got ${ana.S_det}`);
  assert.ok(!ana.branches['c++']?.stable, 'Scenario 2 c++ branch must be unstable by Routh-Hurwitz');

  // Fixed 3.0 us delay: T1 -> T2 locks in m13
  spin.setState({ tone2WaitSS: false, tone2DelayNs: 3.0 });
  spin.fastForwardSteadyState();
  const t1Win13 = spin.cycleAvgModeAmp(13), t1Los16 = spin.cycleAvgModeAmp(16);
  assert.ok(t1Win13 > 0.018 && t1Los16 < 0.001, `T1->T2 (3us) should lock m13 and quench m16: (${t1Win13.toFixed(4)}, ${t1Los16.toFixed(4)})`);

  // Swap T1 <-> T2 (3.0 us delay): T2 -> T1 locks in m16
  spin.swapTonesAction();
  spin.setState({ tone2WaitSS: false, tone2DelayNs: 3.0 });
  spin.fastForwardSteadyState();
  const t2Los13 = spin.cycleAvgModeAmp(13), t2Win16 = spin.cycleAvgModeAmp(16);
  assert.ok(t2Win16 > 0.030 && t2Los13 < 0.001, `Swapped T2->T1 (3us) should lock m16 and quench m13: (${t2Los13.toFixed(4)}, ${t2Win16.toFixed(4)})`);

  // Wait T1 SS automatic release
  spin.applyPreset('nonrecip');
  assert.equal(spin.getState().tone2WaitSS, true);
  spin.fastForwardSteadyState();
  assert.equal(spin.getState().tone2Released, true, 'Wait T1 SS must auto-release Tone 2 after Tone 1 steady state');
  assert.ok(spin.cycleAvgModeAmp(13) > 0.018 && spin.cycleAvgModeAmp(16) < 0.001);
});

test('Q. Experimental AWG (1 ns Tukey, relative phase Δφ), Brown FDT (300 K) & Full NMM DFT', () => {
  spin.setNmmMode(false);

  // 1. 1.0 ns Tukey turn-on ramp: tukey(128, 1)[:64] at 64 GS/s
  assert.equal(spin.TUKEY_RAMP_NS, 1.0);
  spin.setState({ temperatureK: 0, dualTone: true, driveAmp: 0.005, driveAmp2: 0.005, tone2DelayNs: 3.0, tone2WaitSS: false, simTime: 0 });
  assert.ok(Math.abs(spin.tone1EnvAt(0)) < 1e-12);
  assert.ok(Math.abs(spin.tone1EnvAt(spin.nsToTau(0.5)) - 0.5) < 1e-6);
  assert.ok(Math.abs(spin.tone1EnvAt(spin.nsToTau(1.0)) - 1.0) < 1e-6);
  assert.equal(spin.tone2EnvAt(spin.usToTau(2.9)), 0);
  assert.ok(Math.abs(spin.tone2EnvAt(spin.usToTau(3.0)) - 0) < 1e-12);
  assert.ok(Math.abs(spin.tone2EnvAt(spin.usToTau(3.0) + spin.nsToTau(0.5)) - 0.5) < 1e-6);
  assert.ok(Math.abs(spin.tone2EnvAt(spin.usToTau(3.0) + spin.nsToTau(1.0)) - 1.0) < 1e-6);

  // 2. Relative phase dependence (Δφ = 0 deg constructive vs 180 deg destructive on sub-threshold dual tones)
  spin.setState({ damping: 0.002, driveMode: 'parallel' });
  const hSub = 0.75 * spin.parametricThreshold(13);
  spin.setState({
    temperatureK: 0, damping: 0.002, driveMode: 'parallel',
    driveFreq: 2 * spin.modeFreqGHz(13), driveAmp: hSub,
    dualTone: true, driveFreq2: 2 * spin.modeFreqGHz(13), driveAmp2: hSub,
    tone2DelayNs: 0, tone2WaitSS: false, drivePhase1: 0, drivePhase2: 0, simTime: 0
  });
  spin.seedSpins(0.005, [13]);
  spin.fastForwardSteadyState();
  const aInPhase = spin.cycleAvgModeAmp(13);

  spin.setState({ drivePhase2: Math.PI, simTime: 0 });
  spin.seedSpins(0.005, [13]);
  spin.fastForwardSteadyState();
  const aOutPhase = spin.cycleAvgModeAmp(13);
  assert.ok(aInPhase > 0.035 && aOutPhase < 1e-4, `Relative phase must toggle constructive/destructive threshold: 0°=${aInPhase.toFixed(4)} vs 180°=${aOutPhase.toFixed(5)}`);

  // 3. Brown FDT 300 K thermal equipartition & spontaneous nucleation from zero initial state
  spin.setState({ temperatureK: 300, thermalSeed: 12345, damping: 0.002, driveAmp: 0, dualTone: false, driveAmp2: 0, simTime: 0 });
  spin.seedSpins(0.0);
  let sumSq = 0, nSamp = 0;
  for (let i = 0; i < 4000; i++) {
    spin.stepLlg(0.65, true);
    if (i >= 2000 && i % 20 === 0) {
      for (let m = 0; m < spin.M; m++) sumSq += spin.aRe[m] ** 2 + spin.aIm[m] ** 2;
      nSamp += spin.M;
    }
  }
  const rmsSim = Math.sqrt(sumSq / nSamp);
  const rmsTh = spin.thermalRmsAmp();
  assert.ok(rmsSim > 0.8 * rmsTh && rmsSim < 1.25 * rmsTh,
    `FDT 300K equipartition RMS mismatch: sim=${rmsSim.toExponential(3)} vs theory=${rmsTh.toExponential(3)}`);

  // Spontaneous thermal nucleation from exact 0 state under 1.8*h_th
  spin.setState({
    temperatureK: 300, thermalSeed: 777, damping: 0.002, driveMode: 'parallel',
    driveFreq: 2 * spin.modeFreqGHz(13) + 0.002, driveAmp: 1.8 * spin.parametricThreshold(13),
    dualTone: false, simTime: 0
  });
  spin.seedSpins(0.0);
  spin.fastForwardSteadyState();
  assert.ok(spin.modeRmsAmp(13) > 0.04, `300K thermal bath must spontaneously nucleate m13 from a=0: got ${spin.modeRmsAmp(13).toFixed(4)}`);

  // Shot-to-shot stochastic winner selection across thermal seeds in Scenario 2 (13, 16)
  spin.setState({ damping: 0.002, driveMode: 'parallel' });
  const ampStoch13 = 2.4 * spin.parametricThreshold(13);
  const ampStoch16 = 2.5 * spin.parametricThreshold(16);
  let saw13Win = false, saw16Win = false;
  for (const seed of [1, 2]) {
    spin.setState({
      temperatureK: 300, thermalSeed: seed, driveMode: 'parallel', damping: 0.002,
      dualTone: true, tone2WaitSS: false, tone2DelayNs: 0,
      drivePhase1: 0, drivePhase2: 0,
      driveFreq: 2 * spin.modeFreqGHz(13) - 0.001, driveAmp: ampStoch13,
      driveFreq2: 2 * spin.modeFreqGHz(16) - 0.001, driveAmp2: ampStoch16,
      simTime: 0
    });
    spin.seedSpins('thermal');
    spin.fastForwardSteadyState(160000);
    if (spin.modeRmsAmp(13) > 0.05 && spin.modeRmsAmp(16) < 0.005) saw13Win = true;
    if (spin.modeRmsAmp(16) > 0.05 && spin.modeRmsAmp(13) < 0.005) saw16Win = true;
  }
  assert.ok(saw13Win && saw16Win, `Stochastic thermal switching must sample both m13 and m16 attractors across seeds (m13=${saw13Win}, m16=${saw16Win})`);

  // 4. Full NMM (G1+G2+G3+4ABCD) stability + interpolated 0-6 GHz DFT resolution
  spin.setNmmMode(true);
  spin.setState({ temperatureK: 0, damping: 0.002 });
  spin.exciteEigenmode(27, 0.05);
  for (let i = 0; i < 1200; i++) spin.stepLlg(0.65, true);
  spin.updateSpectrum();
  const fsGHz = (2 * Math.PI) / (spin.FREQ_TO_RAD * spin.SPEC_DT);
  let pkK = 1, pkV = 0;
  for (let k = 1; k < spin.SPEC_LEN / 2; k++) {
    if (spin.spectrum[k] > pkV) { pkV = spin.spectrum[k]; pkK = k; }
  }
  const fPeak = (pkK / spin.SPEC_LEN) * fsGHz;
  assert.ok(Math.abs(fPeak - spin.modeFreqGHz(27)) < 0.06,
    `Full NMM DFT peak ${fPeak.toFixed(3)} GHz must match m27=${spin.modeFreqGHz(27).toFixed(3)} GHz`);
  spin.setNmmMode(false);
});

test('R. Slow default playback, dynamic tensor δm = m - m0 profile fidelity & telemetry evaluation logging', () => {
  // 1. Default playback speed is 1× SLOW (0.06 tau/frame) and 0.25× (0.025 tau/frame)
  assert.equal(spin.SPEED_FRAME_STEPS[1][1], 0.06);
  assert.equal(spin.SPEED_FRAME_STEPS[0.25][1], 0.025);
  spin.applyPreset('param');
  assert.equal(spin.getState().playbackSpeed, 1, 'Default playback speed must be gentle 1× SLOW');

  // 2. Dynamic perturbation δm = m - m0 is identically zero at equilibrium (no static rim leakage)
  spin.setState({ temperatureK: 0, driveAmp: 0, dualTone: false, driveAmp2: 0 });
  spin.seedSpins(0.0);
  spin.reconstructSpatialSpins();
  let eqMaxDyn = 0;
  for (let p = 0; p < spin.NACT; p++) {
    const pix = spin.actPix[p];
    const dyn = Math.hypot(spin.my[pix] - spin.m0yFull[pix], spin.mz[pix] - spin.m0zFull[pix]);
    if (dyn > eqMaxDyn) eqMaxDyn = dyn;
  }
  assert.ok(eqMaxDyn < 1e-4, `Equilibrium δm = m - m0 must be ~0 (no 0.384 rim leak), got ${eqMaxDyn}`);

  // 3. Tensor profiles for m13 (bulk standing nodes) vs m27 (uniform Kittel) are distinct & faithful
  const dynProfileY = (mIdx) => {
    spin.exciteEigenmode(mIdx, 0.03);
    const out = new Float64Array(spin.NACT);
    for (let p = 0; p < spin.NACT; p++) out[p] = spin.my[spin.actPix[p]] - spin.m0yFull[spin.actPix[p]];
    return out;
  };
  const p13 = dynProfileY(13), p23 = dynProfileY(23), p27 = dynProfileY(27);
  let mean13 = 0, std13 = 0, mean27 = 0, std27 = 0, dot13_23 = 0, n13 = 0, n23 = 0;
  for (let p = 0; p < spin.NACT; p++) { mean13 += p13[p]; mean27 += p27[p]; }
  mean13 /= spin.NACT; mean27 /= spin.NACT;
  for (let p = 0; p < spin.NACT; p++) {
    std13 += (p13[p] - mean13) ** 2;
    std27 += (p27[p] - mean27) ** 2;
    dot13_23 += p13[p] * p23[p];
    n13 += p13[p] ** 2;
    n23 += p23[p] ** 2;
  }
  const cos13_23 = Math.abs(dot13_23) / Math.sqrt(n13 * n23);
  assert.ok(Math.sqrt(std13 / spin.NACT) > Math.abs(mean13) * 0.5, 'Bulk mode m13 δmy must exhibit strong spatial nodes');
  assert.ok(Math.abs(mean27) > Math.sqrt(std27 / spin.NACT), 'Kittel mode m27 δmy must be spatially quasi-uniform');
  assert.ok(cos13_23 < 0.25, `Tensor modes m13 and m23 δmy profiles must be nearly orthogonal, got cos=${cos13_23.toFixed(3)}`);

  // 4. Moving time-window plot records nonzero mode trajectory and slides smoothly
  spin.applyPreset('coex');
  for (let i = 0; i < 6; i++) {
    spin.stepLlg(0.65, true);
    spin.updatePlotWindow(0.65);
    spin.recordPlotSample();
  }
  spin.drawModeTimePlot();
  const ps = spin.getPlotStats();
  assert.ok(ps.count >= 5 && ps.lastPlotTau > 0 && ps.winTau > 0,
    `Moving time-window plot must record advancing trajectory samples: ${JSON.stringify(ps)}`);

  // 5. Telemetry evaluation logger captures events, trajectory samples & analytical theory comparison
  spin.recordLogEvent('action', 'Test Coex Trigger');
  spin.recordTelemetrySample(true);
  const report = spin.buildEvaluationReport('test_suite');
  assert.equal(report.protocol, 'spin-waves-evaluation-v1');
  assert.equal(report.appVersion, 'spin-waves-nmm-29m-v2');
  assert.ok(report.counts.events >= 2 && report.counts.samples >= 2);
  assert.equal(report.counts.errors, 0);
  assert.equal(report.finalState.theory?.kind, '2mode');
  assert.equal(report.finalState.theory?.scenario, 1);
  assert.ok(Number.isFinite(report.finalState.theory.S_det) && report.finalState.theory.S_det > 0);
});


