/**
 * Companion NMM ODE Verification for spin-waves/reports/idler_math/IDLER_MATH_REPORT.md.
 * Uses scripts/lib/spin-loader.mjs to verify:
 *  1. Tensor provenance: Q2x (inhomogeneous m0), Q2y (transverse 30-deg coupling), P4B, N4A.
 *  2. Sub-threshold 2B/2C non-degenerate idler: RWA Fast (t2B=t2C=false) vs Full NMM (t2B=t2C=true)
 *     under parallel (0 deg) and mixed (transverse / 30 deg) pumping, verifying P_m ~ h1^2 * |a_k|^2.
 *  3. Supra-threshold 4B/4A pre-excitation memory (D > ND) and alpha relaxation scaling.
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { loadSpin } from '../../../scripts/lib/spin-loader.mjs';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const D = loadSpin();
const M = D.M;

function cAbs2(re, im) { return re * re + im * im; }

// 1. Inspect key tensors around modes m17/18 (2.006 GHz), m19 (2.0185), m22 (2.0356), m23 (2.0405), m27 (2.0546), m28 (2.0689)
const modesOfInterest = [17, 18, 19, 22, 23, 27, 28];
const tensorAudit = {
  modes: modesOfInterest.map(m => ({
    m,
    fGHz: +D.fGhzArr[m].toFixed(5),
    gammaGHz_a0002: +D.modeGammaGHz(m, 0.002).toFixed(6),
    hTh_mT_a0002: +D.parametricThresholdMT(m, 0.002).toFixed(4),
    N4A_self: +D.N4A[m * M + m].toFixed(4),
  })),
  offDiagPairs: []
};

for (const h of [17, 19, 22, 23]) {
  for (const k of [23, 27, 28]) {
    if (h === k) continue;
    const idx = (h * M + k) * 2;
    tensorAudit.offDiagPairs.push({
      h, k,
      f_h: +D.fGhzArr[h].toFixed(4),
      f_k: +D.fGhzArr[k].toFixed(4),
      absQ2x: +Math.hypot(D.Q2x[idx], D.Q2x[idx + 1]).toFixed(5),
      absQ2y: +Math.hypot(D.Q2y[idx], D.Q2y[idx + 1]).toFixed(5),
      absP2x: +Math.hypot(D.P2x[idx], D.P2x[idx + 1]).toFixed(5),
      absP2y: +Math.hypot(D.P2y[idx], D.P2y[idx + 1]).toFixed(5),
      absP4B: +Math.hypot(D.P4B[idx], D.P4B[idx + 1]).toFixed(4),
      N4A_hk: +D.N4A[h * M + k].toFixed(4),
    });
  }
}

// 2. Verify sub-threshold 2B/2C idler scaling P_m ~ h1^2 * |a_k|^2 in Full NMM vs RWA-Fast
// Choose target mode k = 28 (2f_28 = 4.1378 GHz) and idler m = 19 (2.0185 GHz) -> f1 = f_19 + f_28 = 4.0874 GHz
function runTwoTone(driveMode, t2BC, h1_mT, h2_mT, f1GHz, f2GHz, alpha = 0.002, t2DelayNs = 0, totalUs = 0.8, seed28 = 0) {
  D.setNmmMode(t2BC ? 'full' : 'rwa');
  D.seedSpins('deterministic');
  D.resetModeHistory();
  if (seed28 > 0) { D.aRe[28] = seed28; D.aIm[28] = 0.0; }
  D.setState({
    temperatureK: 0,
    damping: alpha,
    driveMode,
    driveFreq: f1GHz,
    driveAmp: h1_mT / D.MU0_MS_MT,
    drivePhase1: 0,
    dualTone: true,
    driveFreq2: f2GHz,
    driveAmp2: h2_mT / D.MU0_MS_MT,
    drivePhase2: 0,
    tone2DelayNs: t2DelayNs,
    tone2WaitSS: false,
    terms: { g1: true, t2A: true, t2B: t2BC, t2C: t2BC, g3: false, t4A: true, t4B: true, t4C: false, t4D: false }
  });
  const totalTau = D.usToTau(totalUs);
  const dtTau = 5.0;
  const nSteps = Math.round(totalTau / dtTau);
  const avgStart = Math.round(0.65 * nSteps);
  const acc = new Float64Array(M);
  let cnt = 0;
  for (let s = 0; s < nSteps; s++) {
    D.stepLlg(dtTau);
    if (s >= avgStart) {
      for (let m = 0; m < M; m++) acc[m] += cAbs2(D.aRe[m], D.aIm[m]);
      cnt++;
    }
  }
  for (let m = 0; m < M; m++) acc[m] /= Math.max(1, cnt);
  return { p17: acc[17], p19: acc[19], p23: acc[23], p28: acc[28], pTot: acc.reduce((a, b) => a + b, 0) };
}

const f2Res = 2 * D.fGhzArr[28]; // 4.1378 GHz -> pumps m28
const f1Res = D.fGhzArr[17] + D.fGhzArr[28]; // 4.0751 GHz -> matches m17 + m28 = f1Res with strong Q2x[17,28]=0.07585!
const h2Over = D.parametricThresholdMT(28, 0.002) * 1.9; // strongly supra-threshold Tone 2
const h1Th22 = D.parametricThresholdMT(22, 0.002);

const base0 = runTwoTone('parallel', true, 0, h2Over, f1Res, f2Res, 0.002, 0, 0.8, 0.035);
const idlerSweep = [];
for (const frac of [0.15, 0.30, 0.45, 0.60]) {
  const h1 = frac * h1Th22; // strictly sub-threshold for Tone 1!
  const rRwa = runTwoTone('parallel', false, h1, h2Over, f1Res, f2Res, 0.002, 0, 0.8, 0.035);
  const rPar = runTwoTone('parallel', true, h1, h2Over, f1Res, f2Res, 0.002, 0, 0.8, 0.035);
  const dp17 = rPar.p17 - base0.p17;
  idlerSweep.push({
    h1_mT: +h1.toFixed(4),
    h1_frac: frac,
    rwa_p17: +rRwa.p17.toExponential(4),
    full_p17: +rPar.p17.toExponential(4),
    delta_p17_above_h1_0: +dp17.toExponential(4),
    full_p28: +rPar.p28.toFixed(5),
    delta_p17_over_h1sq_p28: +(dp17 / (h1 * h1 * rPar.p28 + 1e-30)).toExponential(4),
  });
}

// 3. Verify 4B pre-excitation memory (Delayed vs Simultaneous)
const h1Th23 = D.parametricThresholdMT(23, 0.002);
const d_vs_nd = {
  D_delay250ns: runTwoTone('parallel', true, h1Th23 * 1.35, D.parametricThresholdMT(28, 0.002) * 0.96, 2 * D.fGhzArr[23], f2Res, 0.002, 250, 0.8),
  ND_simult:    runTwoTone('parallel', true, h1Th23 * 1.35, D.parametricThresholdMT(28, 0.002) * 0.96, 2 * D.fGhzArr[23], f2Res, 0.002, 0, 0.8),
};

const out = { tensorAudit, idlerSweep, d_vs_nd };
fs.writeFileSync(path.join(__dirname, 'nmm_idler_verification.json'), JSON.stringify(out, null, 2));
console.log(JSON.stringify(out, null, 2));
