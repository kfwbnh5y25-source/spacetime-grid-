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

function loadSpinWavesContext(seedStorage = {}) {
  const noop = () => {};
  const makeCtx = () => new Proxy({
    createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) })
  }, {
    get(target, prop) { return prop in target ? target[prop] : noop; }
  });
  const makeEl = (id = '') => ({
    id,
    style: {},
    classList: { add: noop, remove: noop, toggle: noop },
    setAttribute: noop,
    addEventListener: noop,
    getContext: () => makeCtx(),
    querySelector: () => makeEl(),
    querySelectorAll: () => [],
    closest: () => null,
    scrollIntoView: noop,
    remove: noop,
    contains: () => false,
    dataset: {},
    hidden: id === 'whatsNewModal' || id === 'verifyModal' || id === 'vmTabHandout' || id === 'vmTabLive' || id === 'vmTabDb' || id === 'vmTabExp' || id === 'sweepShell',
    textContent: '',
    innerHTML: '',
    value: 0,
    disabled: false,
    width: 290,
    height: 68
  });
  const elsById = new Map();
  const getEl = id => {
    if (!elsById.has(id)) elsById.set(id, makeEl(id));
    return elsById.get(id);
  };
  const lsStore = new Map(Object.entries(seedStorage).map(([k, v]) => [k, String(v)]));
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
    localStorage: {
      getItem: k => (lsStore.has(k) ? lsStore.get(k) : null),
      setItem: (k, v) => lsStore.set(k, String(v)),
      removeItem: k => lsStore.delete(k)
    },
    fetch: typeof seedStorage.__fetch === 'function' ? seedStorage.__fetch : undefined,
    navigator: { userAgent: 'node-test' },
    setTimeout: () => 0,
    clearTimeout: noop,
    requestAnimationFrame: () => 0,
    addEventListener: noop,
    atob: s => Buffer.from(s, 'base64').toString('binary'),
    document: {
      getElementById: id => getEl(id),
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
    crossPump4BGHz, instantaneousArnoldThresholdMT, computeArnoldHMaxMT, getArnoldRange, arnoldZoomAt, arnoldPanBy, arnoldResetAuto, drawArnoldTongues, phase1Deg, relPhaseDeg,
    cycleAvgModeAmp, modeRmsAmp,
    selfNfsCoeff, mutualNfsCoeff, mutualNfsCoeffGHz, modeNfsShiftGHz, modeMutualNfsShiftGHz,
    steadyStateAmp, analyzeTwoMode, steadyState2Mode, fastForwardSteadyState, fastForwardStepsPerFrame,
    applyPreset, swapTonesAction, SPEED_FRAME_STEPS,
    swLog, recordLogEvent, recordTelemetrySample, buildEvaluationReport,
    RLS_SQL_FIX, getDbConfig, saveDbConfig, resetDbConfig, uploadSwLog, getDbRlsBlockedUrl: () => dbRlsBlockedUrl, formatDbDiagText,
    APP_RELEASE_VERSION, WHATS_NEW_KEY, shouldShowWhatsNew, openWhatsNew, dismissWhatsNew, openVerifyModal, swTutSteps, runSwTutDemo,
    getElHidden: id => !!document.getElementById(id).hidden,
    getElText: id => String(document.getElementById(id).textContent || ''),
    VERIFY_STAGES, runVerificationStage, runFullVerificationProtocol, runThroughputBenchmark,
    inspectLiveCalibration, probe4Avs4BCrossPump,
    G_REF_50PTS, interpolateGrefDB, RF_NOISE_FLOOR, EXP_BCAL_SUMMARY, EXP_ANCHORS_29, EXP_ANCHORS_80_EXTRA, CSI2_ABLATION_DATA,
    awgVpkToUnits, mtToAwgVpk, dualToneRfSummary, syncRfCalUi,
    switchMainTab, applySweepPreset, setSweepMetric, computePointSaSpectrum, computeSingleFreqSweepGrid,
    drawSweepHeatmap, drawSweepSlice, sendSweepPointToLive, exportSweepData,
    getSweepState: () => ({ activeMainTab, ...sweepState }),
    setSweepState: patch => Object.assign(sweepState, patch),
    computeTimeWindowMetrics, buildTimeTraceExport,
    drawModeTimePlot, updatePlotWindow, recordPlotSample,
    getPlotStats: () => ({ count: plotHistCount, winTau: plotWinTau, originTau: plotOriginTau, lastPlotTau, plotYScale, plotYTarget }),
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
      if ('arnoldBand' in patch) arnoldBand = patch.arnoldBand;
      if ('arnoldAutoFollow' in patch) arnoldAutoFollow = patch.arnoldAutoFollow;
      if ('arnoldPanMode' in patch) arnoldPanMode = patch.arnoldPanMode;
      if ('arnoldViewZoom' in patch) arnoldViewZoom = patch.arnoldViewZoom;
      if ('arnoldCenterF' in patch) arnoldCenterF = patch.arnoldCenterF;
      if ('arnoldBaseMT' in patch) arnoldBaseMT = patch.arnoldBaseMT;
      if ('specBand' in patch) specBand = patch.specBand;
      if ('playbackSpeed' in patch) playbackSpeed = patch.playbackSpeed;
      if ('showExpMarkers' in patch) showExpMarkers = patch.showExpMarkers;
      if ('rfChainGainDB' in patch) rfChainGainDB = patch.rfChainGainDB;
      if ('rfKEffMTPerV' in patch) rfKEffMTPerV = patch.rfKEffMTPerV;
      if ('rfVpkTargetMV' in patch) rfVpkTargetMV = patch.rfVpkTargetMV;
      if ('simTime' in patch) { simTime = patch.simTime; resetModeHistory(); recordModeSample(); }
    },
    getState() {
      return {
        driveMode, driveFreq, driveAmp, damping, simTime, dualTone, driveFreq2, driveAmp2,
        tone2DelayNs, tone2WaitSS, tone2Released, tone2ReleaseT,
        drivePhase1, drivePhase2, temperatureK, playbackSpeed, arnoldBand, specBand,
        arnoldAutoFollow, arnoldPanMode, arnoldViewZoom, arnoldCenterF, arnoldBaseMT,
        lastSpecPeakGHz, isRwa: isRwaFast()
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

test('S. Arnold 29-mode distinct thresholds, 4A/4B physics, 300K nonrecip, smooth Y-slew & 1.5–2.3 GHz mode-resolved FFT', () => {
  // 1. All 29 mode thresholds are distinct at alpha=0.002 and unclipped by dynamic ceiling
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, damping: 0.0020, arnoldBand: 'all', specBand: 'modes' });
  let minTh = Infinity, maxTh = 0;
  for (let m = 0; m < spin.M; m++) {
    const th = spin.parametricThresholdMT(m);
    if (th < minTh) minTh = th;
    if (th > maxTh) maxTh = th;
  }
  assert.ok((maxTh / minTh) > 1.25 && (maxTh - minTh) > 0.10, `29-mode thresholds must be distinct (>25% spread), got ${minTh.toFixed(3)}..${maxTh.toFixed(3)} mT`);
  const allRange = spin.getArnoldRange(true);
  assert.ok(Math.abs(allRange.fMin - 3.38) < 1e-6 && Math.abs(allRange.fMax - 4.18) < 1e-6);
  assert.ok(allRange.hMax > maxTh * 1.15, `Arnold ceiling ${allRange.hMax} mT must exceed max tip ${maxTh.toFixed(3)} mT`);
  spin.setState({ damping: 0.0080 });
  let maxTh8 = 0, minTh8 = Infinity;
  for (let m = 0; m < spin.M; m++) {
    const th = spin.parametricThresholdMT(m);
    if (th < minTh8) minTh8 = th;
    if (th > maxTh8) maxTh8 = th;
  }
  const hiRange = spin.getArnoldRange(true);
  assert.ok((maxTh8 - minTh8) > 0.40 && hiRange.hMax > maxTh8 * 1.15, `High-α Arnold ceiling ${hiRange.hMax} mT must unclip max tip ${maxTh8.toFixed(3)} mT`);
  spin.setState({ damping: 0.0020, arnoldBand: 'bulk' });
  const bulkRange = spin.getArnoldRange(true);
  assert.ok(Math.abs(bulkRange.fMin - 3.88) < 1e-6 && Math.abs(bulkRange.fMax - 4.16) < 1e-6);
  spin.setState({ arnoldBand: 'all' });

  // 2. 4A N_hk OFF stops Arnold tongue drift; 4B P_hk lowers threshold & physically cross-pumps m15 from m16
  spin.setState({ simTime: 0 });
  spin.aRe.fill(0); spin.aIm.fill(0); spin.aRe[13] = 0.05;
  spin.terms.t4A = true;
  assert.ok(Math.abs(spin.modeMutualNfsShiftGHz(15)) > 1e-4);
  spin.terms.t4A = false;
  assert.equal(spin.modeNfsShiftGHz(15), 0);
  assert.equal(spin.modeMutualNfsShiftGHz(15), 0);
  spin.terms.t4A = true;

  spin.aRe.fill(0); spin.aIm.fill(0);
  const base15 = spin.instantaneousArnoldThresholdMT(15, 2 * spin.modeFreqGHz(15));
  spin.aRe[16] = 0.045;
  const fTip15 = 2 * (spin.modeFreqGHz(15) + spin.modeMutualNfsShiftGHz(15));
  spin.terms.t4B = false;
  assert.ok(Math.abs(spin.instantaneousArnoldThresholdMT(15, fTip15) - base15) < 1e-6, '4A alone must not alter tip h_th');
  spin.terms.t4B = true;
  const red15 = spin.instantaneousArnoldThresholdMT(15, fTip15);
  assert.ok(red15 < base15 - 0.1, `4B P_hk must lower m15 threshold: base=${base15.toFixed(3)} vs 4B=${red15.toFixed(3)}`);
  const probe4B = spin.probe4Avs4BCrossPump();
  assert.ok(probe4B.thresholdReductionMT > 0.1 && probe4B.m11After4B > 0.004 && probe4B.m11After4AOnly < 0.0016,
    `4B probe failed: ${JSON.stringify(probe4B)}`);
  spin.terms.t4B = false;

  // 3. nonrecip preset at BOTH 0 K and 300 K + relative phase preserved after Swap
  for (const temp of [0, 300]) {
    spin.setState({ temperatureK: temp, thermalSeed: 20260101, drivePhase1: 0, drivePhase2: 0 });
    spin.applyPreset('nonrecip');
    spin.fastForwardSteadyState();
    assert.ok(spin.cycleAvgModeAmp(13) > 0.015 && spin.cycleAvgModeAmp(16) < 0.003,
      `nonrecip T=${temp}K T1->T2 expected m13 win, got (${spin.cycleAvgModeAmp(13).toFixed(4)}, ${spin.cycleAvgModeAmp(16).toFixed(4)})`);
    spin.swapTonesAction();
    assert.equal(spin.relPhaseDeg(), 0);
    spin.fastForwardSteadyState();
    assert.ok(spin.cycleAvgModeAmp(16) > 0.025 && spin.cycleAvgModeAmp(13) < 0.003,
      `nonrecip T=${temp}K Swap expected m16 win, got (${spin.cycleAvgModeAmp(13).toFixed(4)}, ${spin.cycleAvgModeAmp(16).toFixed(4)})`);
  }

  // 4. Non-instantaneous smooth Y-axis rescaling in drawModeTimePlot
  spin.setState({ temperatureK: 0, driveAmp: 0, dualTone: false, driveAmp2: 0, simTime: 0 });
  spin.aRe.fill(0); spin.aIm.fill(0); spin.aRe[13] = 0.12;
  spin.stepLlg(0.5, true);
  spin.recordPlotSample();
  spin.drawModeTimePlot();
  const y1 = spin.getPlotStats();
  assert.ok(y1.plotYScale > 0.0205 && y1.plotYScale < y1.plotYTarget * 0.92,
    `Plot Y scale must slew smoothly instead of snapping: Y1=${y1.plotYScale.toFixed(4)}, target=${y1.plotYTarget.toFixed(4)}`);

  // 5. High-resolution 4096-pt FFT resolves edge mode m0 AND closely-spaced modes (m13=1.9739G vs m16=1.9992G)
  assert.equal(spin.SPEC_LEN, 4096);
  const fsGHz = (2 * Math.PI) / (spin.FREQ_TO_RAD * spin.SPEC_DT);

  // Edge mode m0 visibility at linear amplitude (avoiding 207 GHz Kerr shift at 0.04)
  spin.setState({ damping: 0, temperatureK: 0, driveAmp: 0, dualTone: false, driveAmp2: 0, simTime: 0 });
  spin.exciteEigenmode(0, 0.005);
  spin.setState({ damping: 0, driveAmp: 0, simTime: 0 });
  for (let i = 0; i < 1200; i++) spin.stepLlg(0.65, true);
  spin.updateSpectrum();
  assert.ok(Math.abs(spin.getState().lastSpecPeakGHz - spin.modeFreqGHz(0)) < 0.015,
    `Edge mode m0 FFT peak ${spin.getState().lastSpecPeakGHz} GHz must match ${spin.modeFreqGHz(0).toFixed(4)} GHz`);

  // Dual close modes m13 (1.9739 GHz) + m16 (1.9992 GHz) simultaneously resolved inside 1.5–2.3 GHz
  spin.setState({ damping: 0, temperatureK: 0, driveAmp: 0, dualTone: false, driveAmp2: 0, simTime: 0 });
  spin.seedSpins(0.0);
  spin.aRe[13] = 0.015; spin.aRe[16] = 0.015;
  for (let i = 0; i < 8000; i++) spin.stepLlg(0.65, true);
  spin.updateSpectrum();
  const binAt = fG => Math.round((fG / fsGHz) * spin.SPEC_LEN);
  const k13 = binAt(spin.modeFreqGHz(13)), k16 = binAt(spin.modeFreqGHz(16));
  const peakAround = kCenter => {
    let best = 0;
    for (let k = kCenter - 2; k <= kCenter + 2; k++) if (spin.spectrum[k] > best) best = spin.spectrum[k];
    return best;
  };
  let valley = Infinity;
  for (let k = k13 + 1; k < k16; k++) if (spin.spectrum[k] < valley) valley = spin.spectrum[k];
  const p13 = peakAround(k13), p16 = peakAround(k16);
  assert.ok(p13 > valley * 2.0 && p16 > valley * 2.0,
    `4096-pt FFT must resolve m13 (${spin.modeFreqGHz(13).toFixed(4)}G, pk=${p13.toFixed(2)}) and m16 (${spin.modeFreqGHz(16).toFixed(4)}G, pk=${p16.toFixed(2)}) above valley=${valley.toFixed(2)}`);
});

test('T. Interactive 8-Stage Benchmark, Calibration & Verification Protocol + DB Config', async () => {
  spin.resetDbConfig();
  const cfg0 = spin.getDbConfig();
  assert.ok(cfg0.url.includes('spacetime_feedback_runs'));
  spin.saveDbConfig('https://example.supabase.co/rest/v1/custom', 'jwt.test.token');
  assert.equal(spin.getDbConfig().url, 'https://example.supabase.co/rest/v1/custom');
  spin.resetDbConfig();

  const cal = spin.inspectLiveCalibration();
  assert.ok(cal && Array.isArray(cal.activeModes) && cal.arnoldCeilingMT > 0);

  const ver = await spin.runFullVerificationProtocol();
  assert.equal(ver.stagesCompleted, 8, 'All 8 verification stages must execute');
  assert.ok(ver.totalChecks >= 21, `Expected >=21 protocol checks, got ${ver.totalChecks}`);
  assert.equal(ver.passedChecks, ver.totalChecks,
    `Verification protocol had failing checks: ${JSON.stringify(ver.stages.filter(s => !s.pass), null, 2)}`);

  const rep = spin.buildEvaluationReport('verified_test');
  assert.equal(rep.verification.passedChecks, rep.verification.totalChecks);
  assert.ok(rep.calibration && rep.benchmark && rep.benchmark.rwaKstepsPerSec > 0);
  assert.ok(typeof spin.formatDbDiagText() === 'string');
});

test('U. Unmodified 4B prior-excitation downward Arnold tongue & phase dependence', () => {
  spin.setNmmMode(false);
  spin.terms.t4A = true;
  spin.terms.t4B = true;

  // 1. Coherent phase projection on instantaneousArnoldThresholdMT (with same-tau a_16^2, rotating phi flips wParGHz)
  spin.setState({ temperatureK: 0, damping: 0.0020, simTime: 0 });
  spin.aRe.fill(0); spin.aIm.fill(0);
  spin.aRe[16] = 0.045;
  const f15 = 2 * spin.modeFreqGHz(15) + 0.004;
  const base15 = spin.parametricThresholdMT(15);
  const thOpt = spin.instantaneousArnoldThresholdMT(15, f15, Math.PI / 2);
  const thOpp = spin.instantaneousArnoldThresholdMT(15, f15, -Math.PI / 2);
  assert.ok(thOpt < base15 - 0.15 && thOpp > base15 + 0.15,
    `Phase-resolved 4B threshold must pull DOWN at constructive phase (${thOpt.toFixed(3)} mT) and push UP at opposite phase (${thOpp.toFixed(3)} mT) relative to base=${base15.toFixed(3)} mT`);

  // 2. Sub-threshold cross4b preset (h_p = 0.90 * h_th(15)) at 0 K and 300 K
  for (const temp of [0, 300]) {
    spin.setState({ temperatureK: temp, thermalSeed: 20260101 });
    spin.applyPreset('cross4b');
    assert.equal(spin.phase1Deg(), 0);
    // Without prior m16 excitation, sub-threshold 0.90*h_th(15) decays to floor
    spin.aRe[16] = 0; spin.aIm[16] = 0;
    for (let s = 0; s < 30000; s++) spin.stepLlg(0.5, true);
    const noPrior15 = spin.cycleAvgModeAmp(15);

    // With prior m16 excitation & constructive phase phi_1 = 0 deg, m15 ignites strongly above floor
    spin.setState({ temperatureK: temp, thermalSeed: 20260101 });
    spin.applyPreset('cross4b');
    for (let s = 0; s < 30000; s++) spin.stepLlg(0.5, true);
    const phi0_15 = spin.cycleAvgModeAmp(15);

    // With prior m16 excitation & destructive phase phi_1 = 180 deg, m15 stays suppressed
    spin.setState({ temperatureK: temp, thermalSeed: 20260101 });
    spin.applyPreset('cross4b');
    spin.setState({ drivePhase1: Math.PI });
    assert.equal(spin.phase1Deg(), 180);
    for (let s = 0; s < 30000; s++) spin.stepLlg(0.5, true);
    const phi180_15 = spin.cycleAvgModeAmp(15);

    assert.ok(noPrior15 < 0.0015 && phi0_15 > 0.012 && phi180_15 < 0.003,
      `T=${temp}K cross4b expected sub-threshold ignition only with prior m16 + phi1=0deg: noPrior=${noPrior15.toFixed(4)}, 0deg=${phi0_15.toFixed(4)}, 180deg=${phi180_15.toFixed(4)}`);
    assert.ok(phi0_15 > phi180_15 * 4, `T=${temp}K phase contrast too small: 0deg=${phi0_15.toFixed(4)} vs 180deg=${phi180_15.toFixed(4)}`);
  }
  spin.terms.t4B = false;
});

test('V. Arnold Tongues H/V/Zoom sliders, auto-follow frame lock & trackpad zoom/pan', () => {
  spin.terms.t4A = true;
  spin.setState({ arnoldBand: 'all', driveMode: 'parallel', temperatureK: 0 });
  spin.aRe.fill(0); spin.aIm.fill(0);
  spin.arnoldResetAuto();
  const baseRange = spin.getArnoldRange(true);
  assert.ok(Math.abs(baseRange.fMin - 3.38) < 1e-6 && Math.abs(baseRange.fMax - 4.18) < 1e-6 && baseRange.hMin === 0,
    `Expected default ALL 29 range 3.38..4.18 GHz, got ${baseRange.fMin}..${baseRange.fMax}`);

  // 1. Auto-follow keeps large 4A downward-shifted tongue (e.g. m27 strong softening self-NFS) and high pump crosshair in frame
  spin.aRe[27] = 0.22;
  const f27Shifted = 2 * (spin.modeFreqGHz(27) + spin.modeMutualNfsShiftGHz(27));
  const f27Self = 2 * (spin.modeFreqGHz(27) + spin.modeNfsShiftGHz(27));
  spin.setState({ driveFreq: f27Shifted, driveAmp: 2.2 / spin.MU0_MS_MT });
  const autoRange = spin.getArnoldRange(true);
  assert.ok(f27Shifted >= autoRange.fMin && f27Shifted <= autoRange.fMax,
    `Auto-follow must keep shifted tongue ${f27Shifted.toFixed(4)} GHz in frame [${autoRange.fMin.toFixed(4)}, ${autoRange.fMax.toFixed(4)}]`);
  assert.ok(f27Self < f27Shifted - 0.05 && f27Self >= autoRange.fMin && f27Self <= autoRange.fMax,
    `Auto-follow must also keep single-mode self-NFS tongue ${f27Self.toFixed(4)} GHz in frame [${autoRange.fMin.toFixed(4)}, ${autoRange.fMax.toFixed(4)}]`);
  assert.ok(autoRange.hMax >= 2.2, `Auto-follow hMax (${autoRange.hMax.toFixed(2)}) must cover 2.2 mT drive`);
  spin.drawArnoldTongues();
  assert.ok(spin.getElText('arnoldTermsNote').includes('4A self:'),
    `Expected single-mode self-NFS note in #arnoldTermsNote, got: ${spin.getElText('arnoldTermsNote')}`);

  // 2. Cursor-anchored zoom narrows spans and disables auto-follow
  spin.aRe.fill(0);
  spin.setState({ driveFreq: 3.95, driveAmp: 0.6 / spin.MU0_MS_MT });
  spin.arnoldResetAuto();
  const r0 = spin.getArnoldRange(true);
  const nx = 0.72, ny = 0.25;
  const anchorF = r0.fMin + nx * (r0.fMax - r0.fMin);
  const anchorH = r0.hMin + ny * (r0.hMax - r0.hMin);
  spin.arnoldZoomAt(2.5, nx, ny);
  const rZ = spin.getArnoldRange(true);
  assert.equal(spin.getState().arnoldAutoFollow, false);
  assert.ok(Math.abs((rZ.fMax - rZ.fMin) - (r0.fMax - r0.fMin) / 2.5) < 0.01,
    `Expected frequency span shrunk ~2.5x: ${(rZ.fMax - rZ.fMin).toFixed(4)} vs ${(r0.fMax - r0.fMin).toFixed(4)}`);
  const anchorF2 = rZ.fMin + nx * (rZ.fMax - rZ.fMin);
  const anchorH2 = rZ.hMin + ny * (rZ.hMax - rZ.hMin);
  assert.ok(Math.abs(anchorF2 - anchorF) < 0.02 && Math.abs(anchorH2 - anchorH) < 0.05,
    `Cursor anchor must hold under zoom: F ${anchorF.toFixed(3)}->${anchorF2.toFixed(3)}, H ${anchorH.toFixed(3)}->${anchorH2.toFixed(3)}`);

  // 3. Horizontal & vertical panning move both windows smoothly and stay clamped
  spin.arnoldPanBy(-0.14, 0.35);
  const rP = spin.getArnoldRange(true);
  assert.ok(rP.fMin < rZ.fMin - 0.08 && rP.hMin > rZ.hMin + 0.15,
    `Pan must move both H & V bounds: F ${rZ.fMin.toFixed(3)}->${rP.fMin.toFixed(3)}, hMin ${rZ.hMin.toFixed(2)}->${rP.hMin.toFixed(2)}`);

  // 4. Reset restores AUTO + hMin = 0 + Zoom = 1
  spin.arnoldResetAuto();
  const stR = spin.getState();
  const rReset = spin.getArnoldRange(true);
  assert.ok(stR.arnoldAutoFollow && Math.abs(stR.arnoldViewZoom - 1) < 1e-6 && rReset.hMin === 0,
    `Reset must restore AUTO & Zoom=1: ${JSON.stringify(stR)}`);
});

test('W. FAST->SS smooth frame pacing & high-amplitude adaptive RK4 stability', () => {
  spin.setNmmMode(false); // RWA Fast
  const rwaSteps = spin.fastForwardStepsPerFrame();
  assert.ok(rwaSteps <= 160 && rwaSteps * 0.35 <= 60,
    `RWA FAST->SS per-frame step budget must stay smooth (got ${rwaSteps} steps, ${(rwaSteps * 0.35).toFixed(1)} tau)`);
  spin.setNmmMode(true); // Full NMM
  const fullSteps = spin.fastForwardStepsPerFrame();
  assert.ok(fullSteps <= 64, `Full NMM FAST->SS steps should stay <= 64 (got ${fullSteps})`);

  // High-amplitude stability on m13 (3.5x h_th) in RWA Fast settles cleanly without NaN or mode hopping
  spin.setNmmMode(false);
  const f13 = 2 * spin.modeFreqGHz(13) + 0.004;
  const h13 = 3.5 * spin.parametricThreshold(13);
  spin.setState({ temperatureK: 0, damping: 0.0020, driveMode: 'parallel', dualTone: false, driveFreq: f13, driveAmp: h13, simTime: 0 });
  spin.seedSpins(0.004, [13]);
  spin.fastForwardSteadyState();
  const a13 = spin.cycleAvgModeAmp(13);
  const th13 = spin.steadyStateAmp(13, f13, h13, 0, true);
  assert.ok(Number.isFinite(a13) && Math.abs(a13 - th13) < 0.015,
    `High-amp m13 (3.5*h_th) must remain finite & close to Eq. (33): sim=${a13.toFixed(4)} vs d'Aq=${th13.toFixed(4)}`);

  // Edge mode m0 under strong pump (3.5x h_th) stabilizes smoothly without blowing up
  const f0 = 2 * spin.modeFreqGHz(0) + 0.002;
  const h0 = 3.5 * spin.parametricThreshold(0);
  spin.setState({ temperatureK: 0, damping: 0.0020, driveMode: 'parallel', dualTone: false, driveFreq: f0, driveAmp: h0, simTime: 0 });
  spin.seedSpins(0.004, [0]);
  spin.fastForwardSteadyState();
  const a0 = spin.cycleAvgModeAmp(0);
  assert.ok(Number.isFinite(a0) && a0 > 0.005 && a0 < 0.8,
    `High-amp m0 must remain bounded & stable (got a0=${a0})`);
});

test("X. Long-Term 10-12 us Benchmarks, Stage 8, Handout Tab, SW v79 & Launch What's New Summary", async () => {
  // 1. Stage 8 verification + Handout tab wiring
  assert.equal(spin.VERIFY_STAGES.length, 8);
  assert.equal(spin.VERIFY_STAGES[7].id, 'bench_longterm_steadystate');
  const stage8 = await spin.runVerificationStage(7, false);
  assert.equal(stage8.pass, true, `Stage 8 failed: ${JSON.stringify(stage8)}`);
  assert.equal(stage8.checks.length, 3);

  spin.openVerifyModal('handout');
  assert.equal(spin.getElHidden('verifyModal'), false);
  assert.equal(spin.getElHidden('vmTabHandout'), false);
  assert.equal(spin.getElHidden('vmTabProtocol'), true);
  assert.equal(spin.getElHidden('vmTabLive'), true);
  assert.equal(spin.getElHidden('vmTabExp'), true);
  assert.equal(spin.getElHidden('vmTabDb'), true);

  spin.openVerifyModal('exp');
  assert.equal(spin.getElHidden('vmTabExp'), false);
  assert.equal(spin.getElHidden('vmTabHandout'), true);

  // 2. Launch "What's New" gating (shows once per release version)
  const freshCtx = loadSpinWavesContext();
  assert.equal(freshCtx.APP_RELEASE_VERSION, 'spacetime-shell-v84');
  assert.equal(freshCtx.shouldShowWhatsNew(), true);
  assert.equal(freshCtx.getElHidden('whatsNewModal'), false, 'Fresh launch must auto-open #whatsNewModal');
  freshCtx.dismissWhatsNew(true);
  assert.equal(freshCtx.getElHidden('whatsNewModal'), true);
  assert.equal(freshCtx.shouldShowWhatsNew(), false);

  const seenCtx = loadSpinWavesContext({ [freshCtx.WHATS_NEW_KEY]: freshCtx.APP_RELEASE_VERSION });
  assert.equal(seenCtx.shouldShowWhatsNew(), false);
  assert.equal(seenCtx.getElHidden('whatsNewModal'), true, 'Seen launch must keep #whatsNewModal hidden');
  seenCtx.openWhatsNew();
  assert.equal(seenCtx.getElHidden('whatsNewModal'), false, 'WHAT\'S NEW button must reopen #whatsNewModal');

  // 3. Service Worker v84 + precached BENCHMARK_HANDOUT.md + portrait-safe cedricEggBtn placement
  const swSrc = fs.readFileSync(path.join(rootDir, 'sw.js'), 'utf8');
  assert.ok(swSrc.includes("'spacetime-shell-v84'"), 'sw.js must use spacetime-shell-v84');
  assert.ok(swSrc.includes('./spin-waves/BENCHMARK_HANDOUT.md'), 'sw.js must precache BENCHMARK_HANDOUT.md');
  const htmlSrc = fs.readFileSync(path.join(rootDir, 'spin-waves/index.html'), 'utf8');
  const qbSlice = htmlSrc.slice(htmlSrc.indexOf('<div id="quickLaunchBar">'), htmlSrc.indexOf('<div id="readoutsDock">'));
  assert.ok(qbSlice.includes('id="cedricEggBtn"'), '#cedricEggBtn must live inside #quickLaunchBar so wrapped portrait #topHud never covers it');

  // 4. Saved 10-12 us benchmark JSON & BENCHMARK_HANDOUT.md prints
  const benchJson = JSON.parse(fs.readFileSync(path.join(rootDir, 'benchmarks/longterm-results.json'), 'utf8'));
  assert.equal(benchJson.singleMode10us.length, 29);
  assert.equal(benchJson.pairMatrix10us.length, 406);
  assert.equal(benchJson.nonrecip12us.length, 14);
  assert.equal(benchJson.phaseDual0K.length, 25);
  assert.equal(benchJson.phase4bTransient0K.length, 25);
  assert.equal(benchJson.sustained4bPhase.length, 11);
  assert.equal(benchJson.thresholdCurve.length, 24);
  const p0 = benchJson.phase4bTransient0K.find(r => r.deg === 0);
  const p180 = benchJson.phase4bTransient0K.find(r => Math.abs(r.deg) === 180);
  assert.ok(p0 && p180 && p0.peak15 > p180.peak15 * 8,
    `4B transient phase gate expected >=8x contrast: 0deg=${p0?.peak15}, 180deg=${p180?.peak15}`);
  const th88Off = benchJson.thresholdCurve.find(r => Math.abs(r.mult - 0.88) < 1e-6 && r.t4B === false);
  const th88On = benchJson.thresholdCurve.find(r => Math.abs(r.mult - 0.88) < 1e-6 && r.t4B === true && r.deg === 0);
  assert.ok(th88On && th88Off && th88On.meanTailH > 0.008 && th88Off.meanTailH < 0.0005,
    `10 us continuous 4B expected ON>0.008 vs OFF<0.0005 at 0.88x: ON=${JSON.stringify(th88On)} OFF=${JSON.stringify(th88Off)}`);

  const handoutMd = fs.readFileSync(path.join(rootDir, 'spin-waves/BENCHMARK_HANDOUT.md'), 'utf8');
  assert.ok(handoutMd.includes('3.09911309e+10') && handoutMd.includes('4.7543e-7') && handoutMd.includes('tukey(128,1)[:64]'));
  assert.ok(handoutMd.includes('Verbatim Physical Prints'));
});

test('Y. Experimental RF 50-Ohm Calibration, G_ref Bypass Curve, Mode ID Tables, CSI2 Verdict & Time-Window Stats', () => {
  // 1. 50-Ohm V_pk = V_pp / 2 math & 50-pt G_ref(f) interpolator
  assert.equal(spin.G_REF_50PTS.length, 50);
  const u1V = spin.awgVpkToUnits(1000, 4.11, 0, null);
  assert.equal(u1V.vppMV, 2000);
  assert.ok(Math.abs(u1V.vRmsMV - 1000 / Math.SQRT2) < 1e-9);
  assert.ok(Math.abs(u1V.pAwgMW - 10.0) < 1e-9, `1000 mV_pk into 50 Ohm must equal 10.0 mW, got ${u1V.pAwgMW}`);
  assert.ok(Math.abs(u1V.pAwgDBm - 10.0) < 1e-9, `10 mW must equal +10.0 dBm, got ${u1V.pAwgDBm}`);
  assert.equal(u1V.gRefInBand, true);
  assert.ok(Math.abs(u1V.gRefDB - (-21.33)) < 1e-6, `Expected G_ref(4110 MHz) = -21.33 dB, got ${u1V.gRefDB}`);
  assert.ok(Math.abs(u1V.pSaRefDBm - (-11.33)) < 1e-6);
  assert.equal(u1V.sampleCalibrated, false);
  assert.equal(u1V.hSampleMT, null);

  // Out-of-band clamping check
  const oob = spin.interpolateGrefDB(3.50);
  assert.equal(oob.inBand, false);
  assert.equal(oob.gRefDB, spin.G_REF_50PTS[0]);

  // 2. Uncalibrated vs calibrated k_eff guard (never fabricate sample-plane h_p)
  const uncal = spin.mtToAwgVpk(0.62, null);
  assert.equal(uncal.calibrated, false);
  assert.equal(uncal.vPkMV, null);

  const cal = spin.mtToAwgVpk(0.62, 1.5);
  assert.equal(cal.calibrated, true);
  assert.ok(Math.abs(cal.vPkMV - (0.62 / 1.5) * 1000) < 1e-9);
  const dual = spin.dualToneRfSummary(250, 3.941, 250, 4.110, 0);
  assert.ok(Math.abs(dual.pTotAwgMW - (dual.t1.pAwgMW + dual.t2.pAwgMW)) < 1e-12);
  assert.ok(Math.abs(dual.pTotAwgMW - 1.25) < 1e-9);

  // syncRfCalUi toggles calibrated badge cleanly
  spin.setState({ rfKEffMTPerV: null, rfChainGainDB: 0, rfVpkTargetMV: 250 });
  spin.syncRfCalUi();
  assert.equal(spin.getElHidden('rfSampleWarn'), false);
  spin.setState({ rfKEffMTPerV: 1.5 });
  spin.syncRfCalUi();
  assert.equal(spin.getElHidden('rfSampleWarn'), true);
  spin.setState({ rfKEffMTPerV: null });

  // 3. Linear-power noise floor, 27 mT B_cal retention, 29+80 mode ID & CSI2 verdict
  assert.ok(Math.abs(spin.RF_NOISE_FLOOR.thresh3SigmaDBm - (-70.92)) < 0.02);
  assert.ok(Math.abs(spin.RF_NOISE_FLOOR.threshSnr10DBm - (-66.20)) < 0.02);
  assert.equal(spin.EXP_BCAL_SUMMARY.retainedMT, 27.0);
  assert.ok(Math.abs(spin.EXP_BCAL_SUMMARY.fit29ModesMT - 26.995) < 1e-6);
  assert.equal(spin.EXP_ANCHORS_29.length, 7);
  assert.equal(spin.EXP_ANCHORS_80_EXTRA.length, 3);
  assert.equal(spin.CSI2_ABLATION_DATA.verdict, 'Partially compatible');

  // 4. Time-Dependent Mode Analysis window stats & trace export
  const m13 = 13;
  const f13 = 2 * spin.modeFreqGHz(m13) + 0.003;
  const h13 = 1.6 * spin.parametricThreshold(m13);
  spin.setNmmMode(false);
  spin.setState({ temperatureK: 0, damping: 0.0020, driveMode: 'parallel', dualTone: false, driveFreq: f13, driveAmp: h13, simTime: 0 });
  spin.seedSpins(0.004, [m13]);
  spin.fastForwardSteadyState();
  for (let i = 0; i < 40; i++) {
    for (let s = 0; s < 20; s++) spin.stepLlg(0.35, true);
    spin.updatePlotWindow(20 * 0.35);
    spin.recordPlotSample();
  }
  const ssMetrics = spin.computeTimeWindowMetrics(0.5, 1.0);
  assert.ok(ssMetrics.sampleCount >= 16);
  const row13 = ssMetrics.modes.find(r => r.m === m13);
  assert.ok(row13 && row13.mean > 0.03 && row13.status === 'CONVERGED',
    `Expected steady-state m13 CONVERGED in Time Analysis, got ${JSON.stringify(row13)}`);
  const trace = spin.buildTimeTraceExport();
  assert.equal(trace.model, 'dAquino-2026-NMM-29m');
  assert.equal(trace.B_DC_mT, 27);
  assert.ok(trace.samples.length >= 20 && trace.metrics.sampleCount >= 20);

  // Verify DECAYING detection when pump turns off
  spin.setState({ driveAmp: 0, simTime: 0 });
  spin.seedSpins(0.05, [m13]);
  for (let i = 0; i < 40; i++) {
    for (let s = 0; s < 30; s++) spin.stepLlg(0.35, true);
    spin.updatePlotWindow(30 * 0.35);
    spin.recordPlotSample();
  }
  const decayMetrics = spin.computeTimeWindowMetrics(0, 1);
  const decay13 = decayMetrics.modes.find(r => r.m === m13);
  assert.ok(decay13 && decay13.status === 'DECAYING',
    `Expected unpumped m13 to classify DECAYING, got ${JSON.stringify(decay13)}`);

  // 5. HTML DOM structural integrity: single script block, zero duplicate element IDs
  const htmlSrc = fs.readFileSync(path.join(rootDir, 'spin-waves/index.html'), 'utf8');
  const scriptMatches = htmlSrc.match(/<script\b/g) || [];
  assert.equal(scriptMatches.length, 1, `Expected single <script> tag in spin-waves/index.html, found ${scriptMatches.length}`);
  const idMatches = [...htmlSrc.matchAll(/\bid="([^"]+)"/g)].map(m => m[1]);
  const dupIds = idMatches.filter((id, idx) => idMatches.indexOf(id) !== idx);
  assert.deepEqual(dupIds, [], `Duplicate DOM id(s) found: ${dupIds.join(', ')}`);
  for (const reqId of ['rfCalCard', 'vmTabExp', 'timeEnvBtn', 'timeModeBtn', 'slTimeZoom', 'slTimePan', 'arnoldExpBtn']) {
    assert.ok(idMatches.includes(reqId), `Missing expected UI element #${reqId}`);
  }
});

test('Z. Smooth Dynamic Arnold Tongues at simTime > 0 & Idempotent Supabase RLS 401 Handling', async () => {
  // 1. Dynamic Arnold tongue threshold with fixed driver reference frequency stays smooth at simTime > 0
  spin.setNmmMode(false);
  spin.applyPreset('cross4b');
  const st4b = spin.getState();
  const driveF = st4b.driveFreq;
  const fStart = 3.92, fEnd = 4.02, step = 0.001;
  for (const tauTarget of [0, 2000, 30000]) {
    spin.applyPreset('cross4b');
    while (spin.getState().simTime < tauTarget) spin.stepLlg(0.65, true);
    let maxStepMT = 0;
    let prev = null;
    for (let fp = fStart; fp <= fEnd + 1e-9; fp += step) {
      const hTh = spin.instantaneousArnoldThresholdMT(15, fp, st4b.drivePhase1, driveF);
      assert.ok(Number.isFinite(hTh) && hTh >= 0, `Non-finite threshold at tau=${tauTarget}, fp=${fp}`);
      if (prev !== null) {
        const diff = Math.abs(hTh - prev);
        if (diff > maxStepMT) maxStepMT = diff;
      }
      prev = hTh;
    }
    assert.ok(maxStepMT < 0.25,
      `Dynamic Arnold curve at tau=${tauTarget} must stay smooth (max 1 MHz step=${maxStepMT.toFixed(4)} mT)`);
    assert.doesNotThrow(() => spin.drawArnoldTongues(), `drawArnoldTongues() threw at tau=${tauTarget}`);
  }

  // 2. Idempotent RLS SQL fix + uploadSwLog 401 (42501) resilience & recovery
  assert.ok(spin.RLS_SQL_FIX.includes('DROP POLICY IF EXISTS "allow_anon_insert_runs"'));
  assert.ok(spin.RLS_SQL_FIX.includes('DROP POLICY IF EXISTS "allow_anon_select_runs"'));
  assert.ok(spin.RLS_SQL_FIX.includes('NOTIFY pgrst'));

  let fetchCount = 0;
  let nextStatus = 401;
  const rlsCtx = loadSpinWavesContext({
    __fetch: async () => {
      fetchCount++;
      if (nextStatus === 401) {
        return {
          ok: false,
          status: 401,
          text: async () => JSON.stringify({ code: '42501', message: 'new row violates row-level security policy for table "spacetime_feedback_runs"' })
        };
      }
      return { ok: true, status: 201, text: async () => '' };
    }
  });

  const ok1 = await rlsCtx.uploadSwLog('manual', false);
  assert.equal(ok1, false);
  assert.equal(fetchCount, 1);
  assert.equal(rlsCtx.swLog.dbStatus, 'RLS 401');
  assert.ok(rlsCtx.getDbRlsBlockedUrl() && rlsCtx.getDbRlsBlockedUrl().includes('spacetime_feedback_runs'));

  // Background auto-uploads skip network spam once URL is marked RLS 401
  const okBg = await rlsCtx.uploadSwLog('visibility_hidden', true);
  assert.equal(okBg, false);
  assert.equal(fetchCount, 1, 'Background retry must be suppressed after RLS 401 on same URL');
  assert.equal(rlsCtx.swLog.dbStatus, 'RLS 401');

  // Explicit manual retry still probes the server and clears block when fixed (HTTP 201)
  nextStatus = 201;
  const okFixed = await rlsCtx.uploadSwLog('manual', false);
  assert.equal(okFixed, true);
  assert.equal(fetchCount, 2);
  assert.equal(rlsCtx.swLog.dbStatus, 'DB SAVED ✓');
  assert.equal(rlsCtx.getDbRlsBlockedUrl(), null);

  // 3. Updated 10-step Interactive Guide matches DOM pills, target IDs & executes all demo runners
  assert.equal(spin.swTutSteps.length, 10);
  const htmlSrc = fs.readFileSync(path.join(rootDir, 'spin-waves/index.html'), 'utf8');
  const pills = htmlSrc.match(/class="tutPill"/g) || [];
  assert.equal(pills.length, 10, `Expected 10 .tutPill indicators in DOM, found ${pills.length}`);
  const idSet = new Set([...htmlSrc.matchAll(/\bid="([^"]+)"/g)].map(m => m[1]));
  for (let i = 0; i < spin.swTutSteps.length; i++) {
    const step = spin.swTutSteps[i];
    assert.ok(step.tag.includes(`STEP ${i + 1} OF 10`), `Unexpected tag on step ${i}: ${step.tag}`);
    for (const sel of step.targets || []) {
      const m = sel.match(/^#([A-Za-z0-9_-]+)/);
      if (m) assert.ok(idSet.has(m[1]), `Tutorial step ${i + 1} references missing DOM #${m[1]}`);
    }
    assert.doesNotThrow(() => spin.runSwTutDemo(i), `runSwTutDemo(${i}) threw`);
  }
});

test('AA. Tab 2 Single-Frequency Arnold Sweep Reproduction (Wide, Zoom, SA Metrics, Damping Scaling & Send-to-Tab-1)', () => {
  // 1. Initial state is Tab 1 (live) with #sweepShell hidden; switching to Tab 2 reveals #sweepShell & computes Wide grid
  spin.switchMainTab('live');
  assert.equal(spin.getElHidden('dashShell'), false);
  assert.equal(spin.getElHidden('sweepShell'), true);

  spin.switchMainTab('sweep');
  assert.equal(spin.getElHidden('dashShell'), true);
  assert.equal(spin.getElHidden('sweepShell'), false);
  const st = spin.getSweepState();
  assert.equal(st.activeMainTab, 'sweep');
  assert.equal(st.preset, 'wide');
  assert.ok(st.grid && st.grid.fMinMHz === 3760 && st.grid.fMaxMHz === 4250 && st.grid.vMinMV === 200 && st.grid.vMaxMV === 500);
  assert.ok(st.grid.validFrac > 0.05 && st.grid.validFrac < 0.85,
    `Expected partial excited tongue fraction in Wide sweep, got ${st.grid.validFrac}`);

  // 2. Physical onset: inside m27 Kittel tongue (4110 MHz, 450 mV_pk) strongly exceeds noise floor (200 mV_pk)
  const ptKittel = spin.computePointSaSpectrum(4110, 450, { kEffMTPerV: 1.50, damping: 0.0020, tempK: 300 });
  const ptFloor = spin.computePointSaSpectrum(4110, 200, { kEffMTPerV: 1.50, damping: 0.0020, tempK: 300 });
  const ptBulk13 = spin.computePointSaSpectrum(3948, 450, { kEffMTPerV: 1.50, damping: 0.0020, tempK: 300 });
  assert.ok(ptKittel.signalDB - ptFloor.signalDB > 20,
    `Expected >20 dB parametric Kittel signal jump: peak=${ptKittel.signalDB.toFixed(1)} vs floor=${ptFloor.signalDB.toFixed(1)}`);
  assert.ok(ptBulk13.signalDB - ptFloor.signalDB > 15,
    `Expected >15 dB bulk m13 signal jump: peak=${ptBulk13.signalDB.toFixed(1)} vs floor=${ptFloor.signalDB.toFixed(1)}`);
  assert.equal(ptFloor.thresh10x, false);
  assert.ok(Number.isNaN(ptFloor.widthMHz) && Number.isNaN(ptFloor.nEff), 'Below 10x noise must mask width & n_eff as NaN per plot.py');
  assert.equal(ptKittel.thresh10x, true);
  assert.ok(Number.isFinite(ptKittel.widthMHz) && ptKittel.widthMHz > 0.1 && ptKittel.widthMHz < 3.0);
  assert.ok(Number.isFinite(ptKittel.nEff) && ptKittel.nEff >= 4 && ptKittel.nEff <= 150);

  // 3. Raising Gilbert damping reduces excited tongue fraction & raises threshold
  const gridLowAlpha = spin.computeSingleFreqSweepGrid({ nx: 48, ny: 24, damping: 0.0015, noStore: true });
  const gridHiAlpha = spin.computeSingleFreqSweepGrid({ nx: 48, ny: 24, damping: 0.0045, noStore: true });
  assert.ok(gridHiAlpha.validFrac < gridLowAlpha.validFrac,
    `Higher alpha=0.0045 (${gridHiAlpha.validFrac.toFixed(3)}) must shrink excited tongue vs 0.0015 (${gridLowAlpha.validFrac.toFixed(3)})`);

  // 4. Zoom (4050-4080 MHz) & All29 presets + metric tabs + JSON/CSV export + Send Point -> Tab 1
  spin.applySweepPreset('zoom', true);
  assert.equal(spin.getSweepState().fMinMHz, 4050);
  assert.equal(spin.getSweepState().fMaxMHz, 4080);
  assert.equal(spin.getSweepState().vMinMV, 240);
  assert.equal(spin.getSweepState().vMaxMV, 350);

  spin.setSweepMetric('width');
  assert.doesNotThrow(() => spin.drawSweepHeatmap());
  spin.setSweepMetric('neff');
  assert.doesNotThrow(() => spin.drawSweepHeatmap());
  spin.setSweepMetric('signal');
  assert.doesNotThrow(() => spin.drawSweepSlice());

  const csvOut = spin.exportSweepData('csv');
  assert.ok(csvOut.startsWith('freq_MHz,V_pk_mV,h_pump_mT,signal_dB,width_FWHM_MHz,N_eff,dom_mode'));

  spin.setSweepState({ selF: 4110, selV: 350, kEffMTPerV: 1.50, damping: 0.0020 });
  spin.sendSweepPointToLive();
  assert.equal(spin.getElHidden('dashShell'), false);
  assert.equal(spin.getElHidden('sweepShell'), true);
  const liveSt = spin.getState();
  assert.ok(Math.abs(liveSt.driveFreq - 4.110) < 1e-6);
  assert.ok(Math.abs(liveSt.driveAmp * spin.MU0_MS_MT - 0.350 * 1.50) < 1e-6);
  spin.setState({ rfKEffMTPerV: null });
});


