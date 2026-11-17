// nmm_check.mjs — Drive the NMM solver extracted from spin-waves/index.html from Node.
// All solver functions below marked [VERBATIM] are copied character-for-character
// from index.html (IIFE scope). Globals the functions reference are stubbed at the
// top of this file; every stub is documented. Internal verification only.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const NMM_DATA = JSON.parse(fs.readFileSync(path.join(HERE, 'nmm_data.json'), 'utf8'));

// ---------------------------------------------------------------- data arrays
// [VERBATIM b64ToBuf] index.html ~line 970
function b64ToBuf(b64) {
  const bin = atob(b64);
  const len = bin.length;
  const bytes = new Uint8Array(len);
  for (let i = 0; i < len; i++) bytes[i] = bin.charCodeAt(i);
  return bytes.buffer;
}

// [VERBATIM] index.html lines ~965-970, 981-1000
const M = NMM_DATA.M;
const NX = NMM_DATA.nv, NY = NMM_DATA.nv, N = NX * NY;
const BDC_MT = NMM_DATA.BDC_mT;
const MU0_MS_MT = NMM_DATA.mu0Ms_mT;
const OMEGA_M = NMM_DATA.omegaM;
const F_M = NMM_DATA.fM;
const FREQ_TO_RAD = 1 / F_M;

const fGhzArr = new Float32Array(b64ToBuf(NMM_DATA.fGhz));
const wBarArr = new Float32Array(b64ToBuf(NMM_DATA.wBar));
const norm2Arr = new Float32Array(b64ToBuf(NMM_DATA.norm2));
const qDiagArr = new Float32Array(b64ToBuf(NMM_DATA.qDiag));
const edgeFracArr = new Float32Array(b64ToBuf(NMM_DATA.edgeFrac));
const Kp = new Float32Array(b64ToBuf(NMM_DATA.Kp));
const Km = new Float32Array(b64ToBuf(NMM_DATA.Km));
const bxL = new Float32Array(b64ToBuf(NMM_DATA.bxL));
const byL = new Float32Array(b64ToBuf(NMM_DATA.byL));
const bxA = new Float32Array(b64ToBuf(NMM_DATA.bxA));
const byA = new Float32Array(b64ToBuf(NMM_DATA.byA));

const P2x = new Float32Array(b64ToBuf(NMM_DATA.P2x));
const Q2x = new Float32Array(b64ToBuf(NMM_DATA.Q2x));
const P2y = new Float32Array(b64ToBuf(NMM_DATA.P2y));
const Q2y = new Float32Array(b64ToBuf(NMM_DATA.Q2y));
const N4A = new Float32Array(b64ToBuf(NMM_DATA.N4A));
const P4B = new Float32Array(b64ToBuf(NMM_DATA.P4B));
const C4C = new Float32Array(b64ToBuf(NMM_DATA.C4C));
const D4D = new Float32Array(b64ToBuf(NMM_DATA.D4D));

// -------------------------------------------------------------- stubbed state
// Globals referenced by the copied solver functions. Every stub documented:
// - damping: [VERBATIM value] index.html line ~1092: let damping = 0.0005;
let damping = 0.0005;
// - simTime: [VERBATIM init 0] line ~1100: tau = omega_M * t; we keep 0
let simTime = 0;
// - terms: all toggles forced true per task spec (app default at line 1288 is
//   {g1,t2A,t4A}=true, others false). t2A gates the instantaneous threshold,
//   t4A the NFS shifts, t4B the 4B cross-pump.
const terms = { g1: true, t2A: true, t2B: true, t2C: true, g3: true, t4A: true, t4B: true, t4C: true, t4D: true };
// - ramps: [VERBATIM objects] line ~1106, forced active:false -> ramps inactive
const ramp1 = { active: false, t0: 0, dur: 0, f0: 0, f1: 0, h0: 0, h1: 0, p0: 0, p1: 0 };
const ramp2 = { active: false, t0: 0, dur: 0, f0: 0, f1: 0, h0: 0, h1: 0, p0: 0, p1: 0 };
// - carrier phase anchors: 0 (as at sim start) so channelStateAt returns
//   phi = fDef*FREQ_TO_RAD*tau = 0 at simTime = 0
let carrierAnchorTau1 = 0, carrierAnchorPhi1 = 0, carrierAnchorTau2 = 0, carrierAnchorPhi2 = 0;
let carrierPhase1 = 0, carrierPhase2 = 0;
// - drive state: dualTone true per spec; driveFreq/2 set per experiment below
let driveMode = 'parallel';
let driveFreq = 2 * fGhzArr[13] + 0.002;              // [VERBATIM default] line ~1076
let driveAmp = 0.0034;                                 // [VERBATIM default] line ~1077
let dualTone = true;                                   // forced true per spec (app default false)
let driveFreq2 = 2 * fGhzArr[16];                      // [VERBATIM default] line ~1079
let driveAmp2 = 0.0;                                   // [VERBATIM default] line ~1080
let drivePhase1 = 0.0, drivePhase2 = 0.0;              // [VERBATIM defaults] line ~1086
let tone2DelayNs = 3.0;                                // [VERBATIM default] line ~1081
let tone2WaitSS = true;                                // forced: enables delayed protocol path
// - modal state aRe/aIm: [VERBATIM construction] line ~1114 (Float64Array(M))
const aRe = new Float64Array(M), aIm = new Float64Array(M);
// - buffers used by channelStateAt paths (never on the inactive-ramp branch)
const carrierAnchorPhi1s = 0;
let tone2Released = false, tone2ReleaseT = 0; // only used by tone2EnvAt (not called)

// [VERBATIM] index.html line ~1336 — referenced by crossPump4BGHz; with ramps
// inactive and anchors 0 this branch needs evalRampVal/integrateRampFreq never.
function channelStateAt(ch, tau) {
  const isT2 = ch === 2;
  const r = isT2 ? ramp2 : ramp1;
  const ancTau = isT2 ? carrierAnchorTau2 : carrierAnchorTau1;
  const ancPhi = isT2 ? carrierAnchorPhi2 : carrierAnchorPhi1;
  const fDef = isT2 ? driveFreq2 : driveFreq;
  const hDef = isT2 ? driveAmp2 : driveAmp;
  const pDef = isT2 ? drivePhase2 : drivePhase1;
  const f = r.active ? evalRampVal(r.f0, r.f1, r.t0, r.dur, tau) : fDef;
  const h = r.active ? evalRampVal(r.h0, r.h1, r.t0, r.dur, tau) : hDef;
  const pOff = r.active ? evalRampVal(r.p0, r.p1, r.t0, r.dur, tau) : pDef;
  const phi = (!r.active && ancTau === 0 && ancPhi === 0)
    ? (fDef * FREQ_TO_RAD) * tau
    : ancPhi + FREQ_TO_RAD * integrateRampFreq(ancTau, tau, r, fDef);
  return { freq: f, amp: h, phaseOff: pOff, phi, theta: phi + pOff };
}
// Stubs for the never-evaluated ramp branch of channelStateAt:
function evalRampVal() { throw new Error('evalRampVal called — ramp branch unexpectedly active'); }
function integrateRampFreq() { throw new Error('integrateRampFreq called — ramp branch unexpectedly active'); }

// ------------------------------------------------------- [VERBATIM] functions
// index.html lines ~2168-2199 and ~2250-2400. Copied exactly; only modeGamma's
// neighbours are trimmed where noted.

// returns GHz
function modeFreqGHz(mIdx = 0) { return fGhzArr[mIdx]; }
function modeGammaGHz(mIdx = 0) { return Math.max(1e-7, damping * Kp[(mIdx * M + mIdx) * 2] * F_M); }
function modeGamma(mIdx = 0) { return modeGammaGHz(mIdx) * FREQ_TO_RAD; }
function vParamGHz_mT(mIdx = 0) {
  return 0.5 * Math.hypot(Q2x[(mIdx * M + mIdx) * 2], Q2x[(mIdx * M + mIdx) * 2 + 1]) * F_M / MU0_MS_MT;
}
function parametricThresholdMT(mIdx = 0) {
  return modeGammaGHz(mIdx) / Math.max(1e-9, vParamGHz_mT(mIdx));
}
function parametricThreshold(mIdx = 0) {
  return parametricThresholdMT(mIdx) / MU0_MS_MT;
}
function mutualNfsCoeffGHz(mIdx, jIdx) { return N4A[mIdx * M + jIdx] * F_M; }
function mutualNfsCoeff(mIdx, jIdx) { return mutualNfsCoeffGHz(mIdx, jIdx) * FREQ_TO_RAD; }
function selfNfsCoeff(mIdx) { return mutualNfsCoeff(mIdx, mIdx); }
function modeNfsShiftGHz(mIdx = 0, ignoreToggle = false) {
  if (!ignoreToggle && !terms.t4A) return 0;
  let df = 0;
  for (let j = 0; j < M; j++) {
    const sq = aRe[j] * aRe[j] + aIm[j] * aIm[j];
    if (sq > 0) df += mutualNfsCoeffGHz(mIdx, j) * sq;
  }
  return df;
}
function modeMutualNfsShiftGHz(mIdx = 0, ignoreToggle = false) {
  if (!ignoreToggle && !terms.t4A) return 0;
  let df = 0;
  for (let j = 0; j < M; j++) {
    if (j === mIdx) continue;
    const sq = aRe[j] * aRe[j] + aIm[j] * aIm[j];
    if (sq > 0) df += mutualNfsCoeffGHz(mIdx, j) * sq;
  }
  return df;
}
// Exact 4-magnon parametric cross-pumping rate W_4B = sum_{k != m} P_{mk} a_k^2 from Eq. (9) 4B term
function crossPump4BGHz(mIdx = 0, fProbeGHz = null, pumpPhaseRad = null, pumpRefFreqGHz = null) {
  if (!terms.t4B) return { resRateGHz: 0, totalRateGHz: 0, resRe: 0, resIm: 0, wParGHz: 0, wPerpGHz: 0, cosPhase: 0, maxDriver: -1 };
  let totRe = 0, totIm = 0, resRe = 0, resIm = 0, resMagSum = 0, bestContrib = 0, maxDriver = -1;
  const bw = Math.max(0.008, 3 * modeGammaGHz(mIdx));
  for (let k = 0; k < M; k++) {
    if (k === mIdx) continue;
    const sqk = aRe[k] * aRe[k] + aIm[k] * aIm[k];
    if (sqk < 1e-8) continue;
    const hk2 = (mIdx * M + k) * 2;
    const pR = P4B[hk2] * F_M, pI = P4B[hk2 + 1] * F_M;
    const pAbs = Math.hypot(pR, pI);
    if (pAbs <= 1e-8) continue;
    const k2R = aRe[k] * aRe[k] - aIm[k] * aIm[k];
    const k2I = 2 * aRe[k] * aIm[k];
    const wR = pR * k2R - pI * k2I;
    const wI = pR * k2I + pI * k2R;
    totRe += wR;
    totIm += wI;
    const contrib = pAbs * sqk;
    if (contrib > bestContrib) { bestContrib = contrib; maxDriver = k; }
    if (fProbeGHz !== null) {
      const fkEff = fGhzArr[k] + modeNfsShiftGHz(k);
      const detK = 0.5 * (fProbeGHz - 2 * fkEff);
      const overlap = bw / Math.hypot(bw, detK);
      resRe += wR * overlap;
      resIm += wI * overlap;
      resMagSum += contrib * overlap;
    }
  }
  const totalRateGHz = Math.hypot(totRe, totIm);
  if (fProbeGHz === null) { resRe = totRe; resIm = totIm; resMagSum = totalRateGHz; }
  const resVecMag = Math.hypot(resRe, resIm);
  const resRateGHz = Math.max(resMagSum, resVecMag);
  let wParGHz = resRateGHz, wPerpGHz = 0, cosPhase = 1;
  if (typeof pumpPhaseRad === 'number' && Number.isFinite(pumpPhaseRad) && resVecMag > 1e-9) {
    const i2 = (mIdx * M + mIdx) * 2;
    const qR = Q2x[i2], qI = Q2x[i2 + 1], qMag = Math.hypot(qR, qI);
    if (qMag > 1e-12) {
      const fp = (typeof pumpRefFreqGHz === 'number' && Number.isFinite(pumpRefFreqGHz))
        ? pumpRefFreqGHz
        : (fProbeGHz !== null ? fProbeGHz : 2 * fGhzArr[mIdx]);
      const curF1 = ramp1.active ? channelStateAt(1, simTime).freq : driveFreq;
      const curF2 = ramp2.active ? channelStateAt(2, simTime).freq : driveFreq2;
      const useC1 = !ramp1.active && Math.abs(fp - curF1) < 1e-6;
      const useC2 = dualTone && !ramp2.active && Math.abs(fp - curF2) < 1e-6;
      const thBase = useC1 ? carrierPhase1 : (useC2 ? carrierPhase2 : (fp * FREQ_TO_RAD) * simTime);
      const thP = thBase + pumpPhaseRad;
      // Exact Eq. (9) 2A rotating sideband multiplying a_h*: (-i Q2x / |Q2x|) * e^{i(omega_p tau + phi_p)}
      const c = Math.cos(thP), s = Math.sin(thP);
      const uR = (qI * c + qR * s) / qMag;
      const uI = (qI * s - qR * c) / qMag;
      wParGHz = resRe * uR + resIm * uI;
      wPerpGHz = Math.abs(resRe * uI - resIm * uR);
      cosPhase = wParGHz / resVecMag;
    }
  }
  return { resRateGHz, totalRateGHz, resRe, resIm, wParGHz, wPerpGHz, cosPhase, maxDriver };
}
function arnoldThresholdMT(mIdx, fPumpGHz, fShiftGHz = 0) {
  const hTh = parametricThresholdMT(mIdx);
  const delta = 0.5 * (fPumpGHz - 2 * (fGhzArr[mIdx] + fShiftGHz));
  const gamma = modeGammaGHz(mIdx);
  return hTh * Math.hypot(1, delta / gamma);
}
function instantaneousArnoldThresholdMT(mIdx, fPumpGHz, pumpPhaseRad = null, pumpRefFreqGHz = null) {
  const fMut = modeMutualNfsShiftGHz(mIdx); // strictly 0 if 4A N_hk is unchecked!
  const vm = Math.max(1e-9, vParamGHz_mT(mIdx));
  const gamma = modeGammaGHz(mIdx);
  const delta = 0.5 * (fPumpGHz - 2 * (fGhzArr[mIdx] + fMut));
  const reqRate = Math.hypot(gamma, delta);
  const cp = crossPump4BGHz(mIdx, fPumpGHz, pumpPhaseRad, pumpRefFreqGHz);
  if (typeof pumpPhaseRad === 'number' && Number.isFinite(pumpPhaseRad)) {
    const wVec = Math.hypot(cp.resRe, cp.resIm);
    if (!terms.t2A) return cp.resRateGHz >= reqRate ? 0 : Infinity;
    if (wVec >= reqRate) return 0;
    const radSq = reqRate * reqRate - cp.wPerpGHz * cp.wPerpGHz;
    if (radSq <= 0) return reqRate / vm;
    return Math.max(0, Math.sqrt(radSq) - cp.wParGHz) / vm;
  }
  if (!terms.t2A) return cp.resRateGHz >= reqRate ? 0 : Infinity;
  return Math.max(0, reqRate - cp.resRateGHz) / vm;
}
function arnoldThreshold(mIdx, fPumpGHz, fShiftGHz = 0) {
  return arnoldThresholdMT(mIdx, fPumpGHz, fShiftGHz) / MU0_MS_MT;
}

// d'Aquino (2026) Eq. (33)/(34) single-mode PWM steady-state amplitude
function steadyStateAmp(mIdx, fPumpGHz, hPumpNorm, fMutGHz = 0, pwm = true) {
  const hpMT = hPumpNorm * MU0_MS_MT;
  const gam = modeGammaGHz(mIdx), vm = vParamGHz_mT(mIdx);
  if (hpMT * vm <= gam) return 0;
  const rad = Math.sqrt((vm * hpMT) ** 2 - gam * gam);
  const delta = 0.5 * (fPumpGHz - 2 * (fGhzArr[mIdx] + fMutGHz));
  if (pwm && Math.abs(delta) >= rad) return 0;
  const nkk = mutualNfsCoeffGHz(mIdx, mIdx);
  if (Math.abs(nkk) < 1e-8) return 0;
  const sgn = nkk >= 0 ? 1 : -1;
  if (sgn * delta < -rad) return 0;
  const a2 = (delta + sgn * rad) / nkk;
  return a2 > 0 ? Math.sqrt(a2) : 0;
}

function routhHurwitzStable(c1, c2, c3, c4) {
  if (!(c1 > 0 && c2 > 0 && c3 > 0 && c4 > 0)) return false;
  return (c1 * c2 - c3) > 0 && ((c1 * c2 - c3) * c3 - c1 * c1 * c4) > 0;
}
function checkCoupledStab(m1, m2, a1sq, a2sq, d1, d2, g1, g2) {
  const n11 = mutualNfsCoeffGHz(m1, m1), n12 = mutualNfsCoeffGHz(m1, m2);
  const n21 = mutualNfsCoeffGHz(m2, m1), n22 = mutualNfsCoeffGHz(m2, m2);
  const s1 = n11 * a1sq + n12 * a2sq - d1;
  const s2 = n21 * a1sq + n22 * a2sq - d2;
  const A11 = 4 * s1 * n11 * a1sq, A22 = 4 * s2 * n22 * a2sq;
  const c1 = 2 * (g1 + g2);
  const c2 = 4 * g1 * g2 + A11 + A22;
  const c3 = 2 * g2 * A11 + 2 * g1 * A22;
  const c4 = 16 * s1 * s2 * (n11 * n22 - n12 * n21) * a1sq * a2sq;
  return { stable: routhHurwitzStable(c1, c2, c3, c4), hopf: c4 > 0 && !routhHurwitzStable(c1, c2, c3, c4) };
}

// d'Aquino (2026) Eqs. (55)-(60), Table I Rules 1-12, and Routh-Hurwitz Appendix C.3
function analyzeTwoMode(m1, m2, fp1, hp1Norm, fp2, hp2Norm, pwm = true) {
  const hp1 = hp1Norm * MU0_MS_MT, hp2 = hp2Norm * MU0_MS_MT;
  const n11 = mutualNfsCoeffGHz(m1, m1), n12 = mutualNfsCoeffGHz(m1, m2);
  const n21 = mutualNfsCoeffGHz(m2, m1), n22 = mutualNfsCoeffGHz(m2, m2);
  const g1 = modeGammaGHz(m1), g2 = modeGammaGHz(m2);
  const v1 = vParamGHz_mT(m1), v2 = vParamGHz_mT(m2);
  const above1 = v1 * hp1 > g1, above2 = v2 * hp2 > g2;
  const r1 = above1 ? Math.sqrt((v1 * hp1) ** 2 - g1 * g1) : 0;
  const r2 = above2 ? Math.sqrt((v2 * hp2) ** 2 - g2 * g2) : 0;
  const d1 = 0.5 * (fp1 - 2 * fGhzArr[m1]), d2 = 0.5 * (fp2 - 2 * fGhzArr[m2]);
  const s1p = n11 >= 0 ? 1 : -1, s2p = n22 >= 0 ? 1 : -1;

  const u1Exists = above1 && (s1p * d1 > -r1) && (!pwm || Math.abs(d1) < r1);
  const u1Sq = (above1 && (s1p * d1 > -r1)) ? Math.max(0, (d1 + s1p * r1) / n11) : 0;
  const u1Stab = u1Exists && u1Sq > 0 && (!above2 || Math.abs(d2 - n21 * u1Sq) > r2);

  const u2Exists = above2 && (s2p * d2 > -r2) && (!pwm || Math.abs(d2) < r2);
  const u2Sq = (above2 && (s2p * d2 > -r2)) ? Math.max(0, (d2 + s2p * r2) / n22) : 0;
  const u2Stab = u2Exists && u2Sq > 0 && (!above1 || Math.abs(d1 - n12 * u2Sq) > r1);

  const det = n11 * n22 - n12 * n21;
  const S_det = s1p * s2p * det;
  const branches = {};
  if (above1 && above2 && Math.abs(det) > 1e-12) {
    for (const [lbl, sh, sn] of [['c++', 1, 1], ['c+-', 1, -1], ['c-+', -1, 1], ['c--', -1, -1]]) {
      const b1 = d1 + sh * s1p * r1, b2 = d2 + sn * s2p * r2;
      const a1sq = (b1 * n22 - b2 * n12) / det;
      const a2sq = (n11 * b2 - n21 * b1) / det;
      if (a1sq > 0 && a2sq > 0) {
        const st = checkCoupledStab(m1, m2, a1sq, a2sq, d1, d2, g1, g2);
        branches[lbl] = { a1: Math.sqrt(a1sq), a2: Math.sqrt(a2sq), stable: st.stable, hopf: st.hopf };
      }
    }
  }

  // Select predicted outcome under configured Tone 1 -> Tone 2 protocol (or simultaneous)
  const isDelayed = tone2WaitSS || tone2DelayNs > 0;
  let pred = [0, 0], regime = 'z';
  if (isDelayed) {
    if (Math.abs(d1) < r1 && above1 && u1Sq > 0) {
      const a1First = Math.sqrt(u1Sq);
      if (!above2 || Math.abs(d2 - n21 * u1Sq) >= r2) {
        pred = [a1First, 0]; regime = u1Stab ? 'uh' : 'uh*';
      } else if (S_det > 0 && branches['c++'] && branches['c++'].stable) {
        pred = [branches['c++'].a1, branches['c++'].a2]; regime = 'c++';
      } else if (S_det < 0 && branches['c+-'] && branches['c+-'].stable) {
        pred = [branches['c+-'].a1, branches['c+-'].a2]; regime = 'c+-';
      } else if (u2Stab && !u1Stab) {
        pred = [0, Math.sqrt(u2Sq)]; regime = 'un';
      } else {
        pred = [a1First, 0]; regime = 'Q/osc';
      }
    } else if (Math.abs(d2) < r2 && above2 && u2Sq > 0) {
      pred = [0, Math.sqrt(u2Sq)]; regime = 'un';
    }
  } else {
    if (branches['c++'] && branches['c++'].stable) {
      pred = [branches['c++'].a1, branches['c++'].a2]; regime = 'c++';
    } else if (S_det < 0 && branches['c+-'] && branches['c+-'].stable) {
      pred = [branches['c+-'].a1, branches['c+-'].a2]; regime = 'c+-';
    } else if (S_det < 0 && branches['c-+'] && branches['c-+'].stable) {
      pred = [branches['c-+'].a1, branches['c-+'].a2]; regime = 'c-+';
    } else if (u1Stab && !u2Stab) {
      pred = [Math.sqrt(u1Sq), 0]; regime = 'uh';
    } else if (u2Stab && !u1Stab) {
      pred = [0, Math.sqrt(u2Sq)]; regime = 'un';
    } else if (u1Stab && u2Stab) {
      const z1 = Math.sqrt(Math.max(0, (v1 * hp1) ** 2 - d1 * d1)) - g1;
      const z2 = Math.sqrt(Math.max(0, (v2 * hp2) ** 2 - d2 * d2)) - g2;
      if (z2 > z1) { pred = [0, Math.sqrt(u2Sq)]; regime = 'uh|un(n)'; }
      else { pred = [Math.sqrt(u1Sq), 0]; regime = 'uh|un(h)'; }
    }
  }
  return {
    pred, regime, S_det, scenario: S_det >= 0 ? 1 : 2,
    u1: u1Sq > 0 ? Math.sqrt(u1Sq) : 0, u1Stab,
    u2: u2Sq > 0 ? Math.sqrt(u2Sq) : 0, u2Stab,
    branches
  };
}

function steadyState2Mode(i, j, fp1, hp1Norm, fp2, hp2Norm) {
  return analyzeTwoMode(i, j, fp1, hp1Norm, fp2, hp2Norm, true).pred;
}

function closestArnoldMode(fPumpGHz) {
  let bestM = 13, bestD = 1e9;
  for (let m = 0; m < M; m++) {
    const d = Math.abs(fPumpGHz - 2 * fGhzArr[m]);
    if (d < bestD) { bestD = d; bestM = m; }
  }
  return bestM;
}

// ================================================================= ANALYSIS
const out = [];
function P(s) { out.push(s); }
const fmt = (x, d = 6) => (Number.isFinite(x) ? x.toFixed(d) : x.toString());

// Partial verify of verbatim copy against source lines (cheap guard):
P('# NMM solver cross-check — machine-readable stdout');
P('# damping = ' + damping + ' (verbatim index.html line 1092), terms all true, simTime=0, dualTone=true, tone2WaitSS=true');
P('# M=' + M + ', F_M=' + F_M + ' GHz, MU0_MS_MT=' + MU0_MS_MT + ' mT, BDC=' + BDC_MT + ' mT');

// ---- 1. Mode table
P('\n== [1] MODE TABLE (29 modes) ==');
P('mIdx\tf_m[GHz]\tGamma_m[GHz]\tGamma_m[MHz]\tv_m[GHz/mT]\thTh0[mT]\tf^2f_ang\t2f_m[GHz]');
const gamArr = [], vArr = [];
for (let m = 0; m < M; m++) {
  const g = modeGammaGHz(m), v = vParamGHz_mT(m);
  gamArr.push(g); vArr.push(v);
  P([m, fGhzArr[m].toFixed(6), g.toFixed(9), (g * 1000).toFixed(4), v.toFixed(6),
     parametricThresholdMT(m).toFixed(6), Math.abs(fGhzArr[m] - F_M / 2).toFixed(6),
     (2 * fGhzArr[m]).toFixed(6)].join('\t'));
}
P('Gamma range (GHz): ' + Math.min(...gamArr).toFixed(9) + ' .. ' + Math.max(...gamArr).toFixed(9));
P('Gamma range (MHz): ' + (Math.min(...gamArr) * 1000).toFixed(4) + ' .. ' + (Math.max(...gamArr) * 1000).toFixed(4));
P('NOTE: solver Gamma = damping * Kp_diag * F_M, damping = 0.0005. Doc quotes 8.2-9.8 MHz at alpha=0.002 (4x larger).');

// ---- 2. Single-mode Arnold threshold vs Eq.(2)
P('\n== [2] SINGLE-MODE ARNOLD THRESHOLD: arnoldThresholdMT vs Eq.(2) h_th = Gamma*hypot(1, delta/Gamma)/V ==');
P('(Gamma*hypot(1, delta/Gamma) alone is a RATE in GHz; the threshold in mT needs the 1/V. Both compared below.)');
const reps = [];
for (let m = 0; m < M; m++) if (fGhzArr[m] >= 1.9 && fGhzArr[m] <= 2.1) reps.push(m);
P('Modes with f_m in [1.9, 2.1] GHz: ' + reps.join(','));
const detunings = [0, 3, 6, 10]; // MHz
let maxRelDev = 0, maxRelDevInfo = '', maxRelDevRate = 0;
P('mIdx\tf_m[GHz]\tdelta[MHz]\tfPump[GHz]\tarnoldThresholdMT[mT]\tEq2[mT]\trel.devField\tv*arnold[GHz]\tEq2rate[GHz]\trel.devRate');
for (const m of reps) {
  const g = modeGammaGHz(m);
  for (const dMHz of detunings) {
    const fp = 2 * fGhzArr[m] + 2 * dMHz / 1000;
    const solver = arnoldThresholdMT(m, fp);
    const eq2Rate = g * Math.hypot(1, dMHz / 1000 / g); // sqrt(G^2+d^2) in GHz (doc Eq.2 numerator)
    const eq2 = eq2Rate / vArr[m];                      // Eq.(2) full, mT
    const rel = Math.abs(solver - eq2) / eq2;
    const relR = Math.abs(vArr[m] * solver - eq2Rate) / eq2Rate;
    if (rel > maxRelDev) { maxRelDev = rel; maxRelDevInfo = `m=${m} delta=${dMHz}`; }
    if (relR > maxRelDevRate) maxRelDevRate = relR;
    P([m, fGhzArr[m].toFixed(6), dMHz, fp.toFixed(7), solver.toFixed(9), eq2.toFixed(9), rel.toExponential(2),
       (vArr[m] * solver).toExponential(6), eq2Rate.toExponential(6), relR.toExponential(2)].join('\t'));
  }
}
P('MAX RELATIVE DEVIATION (field units, solver vs full Eq.(2)) = ' + maxRelDev.toExponential(3) + ' at ' + maxRelDevInfo);
P('MAX RELATIVE DEVIATION (rate units, v*arnoldThresholdMT vs Gamma*hypot(1,delta/Gamma)) = ' + maxRelDevRate.toExponential(3));
P('Also solver hTh0 vs Gamma/v: ' + [reps[0]].map(m => {
  const mt = parametricThresholdMT(m), gamV = modeGammaGHz(m) / vParamGHz_mT(m);
  return `m=${m}: parametricThresholdMT=${mt.toFixed(9)} mT, Gamma/v=${gamV.toFixed(9)} mT, diff=${(mt - gamV).toExponential(2)}`;
}).join('; '));

// ---- 3. Phase-dependent two-mode threshold at experimental pumps
P('\n== [3] PHASE-DEPENDENT TWO-MODE THRESHOLD (experimental pumps) ==');
const FP1 = 4.070, FP2 = 4.142; // experimental pump frequencies (GHz)
const m1 = closestArnoldMode(FP1), m2 = closestArnoldMode(FP2);
const g1 = modeGammaGHz(m1), g2 = modeGammaGHz(m2);
const v1 = vParamGHz_mT(m1), v2 = vParamGHz_mT(m2);
P(`closestArnoldMode(${FP1}) = mode idx ${m1} (solution idx ${NMM_DATA.selIdx[m1]}), f_m = ${fGhzArr[m1].toFixed(6)} GHz, match |fp-2fm| = ${Math.abs(FP1 - 2 * fGhzArr[m1]).toFixed(6)} GHz`);
P(`closestArnoldMode(${FP2}) = mode idx ${m2} (solution idx ${NMM_DATA.selIdx[m2]}), f_m = ${fGhzArr[m2].toFixed(6)} GHz, match |fp-2fm| = ${Math.abs(FP2 - 2 * fGhzArr[m2]).toFixed(6)} GHz`);
P(`Gamma1 = ${(g1 * 1000).toFixed(4)} MHz, Gamma2 = ${(g2 * 1000).toFixed(4)} MHz; v1 = ${v1.toFixed(6)}, v2 = ${v2.toFixed(6)} GHz/mT`);

// pump1 amplitude such that v1*h1 = 1.3*Gamma1 : hp1_norm = 1.3*g1/(v1*MU0_MS_MT)
const hp1Norm = 1.3 * g1 / (v1 * MU0_MS_MT);
P(`hp1 fixed: v1*h1 = 1.3*Gamma1 -> h1 = ${fmt(1.3 * g1 / v1)} mT = ${hp1Norm.toExponential(6)} (normalized by MU0_MS_MT)`);
// steady state of seed mode (tone 1 below ITS threshold by construction: v*h = 1.3 G > G -> actually above 1-photon threshold!)
// NOTE: v*h1 = 1.3*Gamma1 > Gamma1 means pump1 ALONE is above single-mode threshold.
// The doc scenario needs seed BELOW parametric ignition; here we follow task spec literally.
const a1 = steadyStateAmp(m1, FP1, hp1Norm, 0, true);
P(`steadyStateAmp(seed, FP1, hp1Norm, fMut=0, pwm=true) = ${a1.toExponential(6)}  (|a|^2 = ${(a1 * a1).toExponential(6)})`);
const a1Above = v1 * hp1Norm * MU0_MS_MT > g1;
P(`seed single-mode condition: v1*h1 = ${(v1 * hp1Norm * MU0_MS_MT).toExponential(4)} GHz vs Gamma1 = ${(g1).toExponential(4)} GHz -> above threshold: ${a1Above}`);

// Seed ONLY the seed mode in aRe/aIm for the cross-pump evaluation
aRe.fill(0); aIm.fill(0);
aRe[m1] = a1;

const Wcp = crossPump4BGHz(m2, FP2, null, FP2); // no phase -> total
P(`crossPump4BGHz(target, FP2) (no phase): totalRateGHz = ${Wcp.totalRateGHz.toExponential(6)} GHz, resRateGHz = ${Wcp.resRateGHz.toExponential(6)} GHz, resRe = ${Wcp.resRe.toExponential(6)}, resIm = ${Wcp.resIm.toExponential(6)}, maxDriver idx = ${Wcp.maxDriver} (|a_k|^2 = ${(aRe[Wcp.maxDriver] ** 2).toExponential(3)})`);
const W = Math.hypot(Wcp.resRe, Wcp.resIm);
P(`W (rate) = ${W.toExponential(6)} GHz; W/v2 = ${(W / v2).toExponential(6)} mT; W/Gamma2 = ${(W / g2).toFixed(6)}`);
const chi = Math.atan2(Wcp.resIm, Wcp.resRe); // absolute phase of internal pump vector in solver convention
const q2i2 = (m2 * M + m2) * 2, q2r = Q2x[q2i2], q2i = Q2x[q2i2 + 1];
P(`internal pump angle chi = atan2(resIm,resRe) = ${chi.toFixed(6)} rad; arg(Q2x_diag(target)) = ${Math.atan2(q2i, q2r).toFixed(6)} rad`);

// NFS shift seen by target because seed is excited:
const fMut2 = modeMutualNfsShiftGHz(m2);
const fShift2 = modeNfsShiftGHz(m2);
P(`modeMutualNfsShiftGHz(target) = ${fMut2.toExponential(6)} GHz (self N4A excluded), modeNfsShiftGHz(target) = ${fShift2.toExponential(6)} GHz (includes self)`);
const selfNkk = mutualNfsCoeffGHz(m2, m2);
P(`selfNfsCoeffGHz(target) = ${(selfNkk).toExponential(6)} GHz per |a|^2`);
const d2_solv = 0.5 * (FP2 - 2 * (fGhzArr[m2] + fMut2)); // what instantaneousArnoldThresholdMT uses
const d2_doc = 0.5 * (FP2 - 2 * fGhzArr[m2]);
P(`delta2 used by solver (incl. mutual NFS): ${d2_solv.toExponential(6)} GHz (${(d2_solv * 1000).toFixed(4)} MHz); doc delta2 (bare): ${d2_doc.toExponential(6)} GHz (${(d2_doc * 1000).toFixed(4)} MHz)`);

// 16 relative phases of pump2, fixed hp2 (choose below bare threshold at delta2):
// take hp2 such that v2*h2 = 0.6*reqRate0 where reqRate0 = hypot(Gamma2, d2_bare)
const req0 = Math.hypot(g2, Math.abs(d2_doc));
const hp2Norm = 0.6 * req0 / (v2 * MU0_MS_MT);
P(`hp2 fixed: v2*h2 = 0.6*req0 = ${fmt(0.6 * req0)} GHz -> h2 = ${fmt(0.6 * req0 / v2)} mT = ${hp2Norm.toExponential(6)} normalized`);

const NP = 16;
P('\n-- threshold vs pump2 phase (instantaneousArnoldThresholdMT at FP2, ref FP2) --');
P('psiIdx\tpsi[rad]\thTh_solver[mT]\tWpar[GHz]\tWperp[GHz]');
const solv = [];
for (let i = 0; i < NP; i++) {
  const psi = 2 * Math.PI * i / NP;
  const h = instantaneousArnoldThresholdMT(m2, FP2, psi, FP2);
  const cp = crossPump4BGHz(m2, FP2, psi, FP2);
  solv.push({ psi, h, wPar: cp.wParGHz, wPerp: cp.wPerpGHz });
  P([i, psi.toFixed(6), Number.isFinite(h) ? h.toFixed(9) : 'Infinity', cp.wParGHz.toExponential(6), cp.wPerpGHz.toExponential(6)].join('\t'));
}

// Doc Eq.(6): h_th(Psi) = [sqrt(G2^2 + d2^2 - W^2 sin^2(Psi)) - W cos(Psi)]/v2 — this IS the
// exact solver formula (task labels it Eq.(7); Eq.(7) in the doc is its linearization).
// The only freedom is a constant phase offset (convention of the internal pump angle).
// Fit offset phi0 minimizing max |h_solver(psi) - h_doc(psi - phi0)|, then report residuals.
function docEq6(psi, Wv, d, G, v) {
  const s2 = Wv * Wv * Math.sin(psi) ** 2;
  const rad = G * G + d * d - s2;
  const r = rad < 0 ? 0 : Math.sqrt(rad);
  return Math.max(0, (r - Wv * Math.cos(psi))) / v;
}
function docEq7lin(psi, Wv, d, G, v) {
  return (Math.hypot(G, d) - Wv * Math.cos(psi)) / v;
}
function bestOffset(fn, solvArr, Wv, d, G, v) {
  let best = { phi0: 0, maxDev: Infinity, rms: Infinity };
  for (let k = 0; k < 720; k++) {
    const phi0 = 2 * Math.PI * k / 720;
    let maxDev = 0, ssum = 0, n = 0;
    for (const s of solvArr) {
      if (!Number.isFinite(s.h)) continue;
      const hd = fn(s.psi - phi0, Wv, d, G, v);
      const dev = Math.abs(s.h - hd);
      maxDev = Math.max(maxDev, dev); ssum += dev * dev; n++;
    }
    const rms = Math.sqrt(ssum / Math.max(1, n));
    if (maxDev < best.maxDev) best = { phi0, maxDev, rms };
  }
  return best;
}
// d used by solver includes mutual NFS shift of the target from the seed:
{
  const Wv = W, d = d2_solv, G = g2, v = v2;
  const b6 = bestOffset(docEq6, solv, Wv, d, G, v);
  const b7 = bestOffset(docEq7lin, solv, Wv, d, G, v);
  P(`\nDoc Eq.(6) exact fit (W=${Wv.toExponential(4)} GHz, d=${d.toExponential(4)} GHz, G=${G.toExponential(4)} GHz): fitted internal-phase offset phi0 = ${b6.phi0.toFixed(4)} rad, max|h_solver-h_doc| = ${b6.maxDev.toExponential(3)} mT, rms = ${b6.rms.toExponential(3)} mT over ${NP} phases`);
  P(`Doc Eq.(7) linearized fit (same W,d,G): phi0 = ${b7.phi0.toFixed(4)} rad, max dev = ${b7.maxDev.toExponential(3)} mT, rms = ${b7.rms.toExponential(3)} mT`);
  P(`peak-to-valley of solver curve: ${(Math.max(...solv.map(s => s.h)) - Math.min(...solv.filter(s => Number.isFinite(s.h)).map(s => s.h))).toFixed(6)} mT; 2W/v2 = ${fmt(2 * Wv / v)} mT`);
  P(`chi - phi0 relationship: phase of internal pump predicted = chi - phi0 = ${(chi - b6.phi0).toFixed(4)} rad (doc Psi = phi - phi_W)`);
}
// Also with d2 bare (no NFS) to show what a mismatch means:
{
  const b6b = bestOffset(docEq6, solv, W, d2_doc, g2, v2);
  P(`Doc Eq.(6) with BARE d2 (ignores NFS from seed): best max dev = ${b6b.maxDev.toExponential(3)} mT (larger => NFS matters)`);
}
// phase-averaged check
{
  const avg = solv.reduce((s, x) => s + x.h, 0) / NP;
  const b0 = Math.hypot(g2, d2_solv) / v2;
  P(`phase-averaged solver threshold = ${avg.toFixed(6)} mT vs bare (2A only) = ${b0.toFixed(6)} mT; doc: unchanged to linear order.`)
}

// analyzeTwoMode regime across 16 phases (delayed path), fixed hp1/hp2, experimental pumps
P('\n-- analyzeTwoMode regime vs pump2 phase (delayed protocol: tone2WaitSS=true) --');
P('psiIdx\tpsi[rad]\tregime\tpred[a1,a2]');
for (let i = 0; i < NP; i++) {
  const psi = 2 * Math.PI * i / NP;
  // note: analyzeTwoMode has no phase knob; phase enters only via crossPump4B. Report regime grid anyway.
  const a = analyzeTwoMode(m1, m2, FP1, hp1Norm, FP2, hp2Norm, true);
  P([i, psi.toFixed(6), a.regime, a.pred.map(x => x.toExponential(4)).join(',')].join('\t'));
}

// ---- 4. Delayed-tone protocol: regimes vs pump2 detuning AND pump2 amplitude
P('\n== [4] DELAYED-TONE PROTOCOL: regime vs (Delta, hp2) (tone2WaitSS=true, seed at steady state) ==');
P(`pump1 fixed at FP1 = ${FP1} GHz (seed idx ${m1}), hp1 = 1.3*G1/v1 (normalized ${hp1Norm.toExponential(4)}), |a_seed|^2 = ${(a1 * a1).toExponential(6)}`);
P(`target engaged in delayed path iff |d2 - n21*|a1|^2| < r2 with n21*|a1|^2 = ${(mutualNfsCoeffGHz(m2, m1) * a1 * a1 * 1000).toFixed(4)} MHz (NFS shift of target by seed)`);
P('n21 (GHz per |a|^2) = ' + mutualNfsCoeffGHz(m2, m1).toExponential(6) + '; d2_bare grows linearly in Delta: d2 = (FP1+Delta)/2 - f_target');
P('closestArnoldMode switch points (seed idx fixed ' + m1 + ') in Delta: between 2f targets -> fp2 midpoints:');
for (let t = 0; t < M - 1; t++) {
  const mid = (2 * fGhzArr[t] + 2 * fGhzArr[t + 1]) / 2;
  if (mid > FP1 + 0.030 && mid < FP1 + 0.110) P(`  fp2 = ${mid.toFixed(6)} GHz -> Delta = ${((mid - FP1) * 1000).toFixed(2)} MHz (target idx ${t} -> ${t + 1})`);
}
P('\nREGIME MAP rows = Delta [MHz], cols = hp2 [mT] {0.10, 0.14, 0.20, 0.30, 0.50, 0.80, 1.50}:');
const h2List = [0.10, 0.14, 0.20, 0.30, 0.50, 0.80, 1.50];
P('Delta[MHz]\ttgt\t|' + h2List.slice(1).map(() => ' ').join('') + 'regimes for hp2 = 0.10 .. 1.50 mT|  (delayed protocol)');
for (let DM = 40; DM <= 100; DM += 4) {
  const fp2x = FP1 + DM / 1000;
  const tx = closestArnoldMode(fp2x);
  const regimes = h2List.map(h2mT => {
    tone2WaitSS = true; tone2DelayNs = 0;
    const a = analyzeTwoMode(m1, tx, FP1, hp1Norm, fp2x, h2mT / MU0_MS_MT, true);
    tone2WaitSS = true; tone2DelayNs = 3.0;
    return a.regime;
  });
  P([DM, tx, ...regimes].join('\t'));
}
P('\nBoundary column: distance of d2 from the lock window edge at hp2 = 0.14 mT (target idx ' + m2 + '):');
P('Delta[MHz]\td2bare[MHz]\tn21|a1|^2[MHz]\tr2(hp2=0.14mT)[MHz]\t|d2-n21a1^2|-r2 [MHz]  (>0 => regime uh, seed-only)');
const hp2ref = 0.14 / MU0_MS_MT;
for (let DM = 40; DM <= 100; DM += 4) {
  const fp2x = FP1 + DM / 1000;
  const tx = closestArnoldMode(fp2x);
  const dx = 0.5 * (fp2x - 2 * fGhzArr[tx]);
  const n21a1 = mutualNfsCoeffGHz(tx, m1) * a1 * a1;
  const r2c = Math.sqrt(Math.max(0, (vParamGHz_mT(tx) * hp2ref * MU0_MS_MT) ** 2 - modeGammaGHz(tx) ** 2));
  P([DM, (dx * 1000).toFixed(4), (n21a1 * 1000).toFixed(4), (r2c * 1000).toFixed(4), ((Math.abs(dx - n21a1) - r2c) * 1000).toFixed(3)].join('\t'));
}

// memory-lifetime relevance probe: time scale the solver does NOT own.
P('\n-- What time scales does the solver contain? --');
P('The delayed path is a STATIC SS-selection rule: it never integrates W(t_gap)=W0*exp(-2*Gamma_k*t_gap).');
P('Omega_1 (seed ringdown) = 1/Gamma1 = ' + fmt(1 / g1, 3) + ' ns-equivalent units: Gamma1=' + (g1 * 1000).toFixed(4) + ' MHz -> 1/Gamma1 = ' + fmt(1000 / (g1 * 1000), 4) + ' ns; tau_mem=1/(2 Gamma1) = ' + fmt(500 / (g1 * 1000), 4) + ' ns');
P('i.e. with damping=0.0005 the doc Eq.(10) predicts tau_mem ~ ' + fmt(500 / (g1 * 1000), 3) + ' ns, NOT a function of Delta: no regime boundary in analyzeTwoMode can produce a Delta-dependent memory lifetime');

// convention checks: possible factor-2 spots
P('\n== [5] CONVENTION / FACTOR CHECKS ==');
// (a) vParam: solver v = 0.5*|Q2x_diag| * F_M / MU0_MS_MT  [GHz/mT]
//      growth matrix element = v*h ; doc V defined so off-diag pump element = i V h. Same side.
// (b) delta: solver delta = 0.5*(fp - 2 fm)  == doc delta = 2pi(f_m - fp/2)/(2pi) sign-flipped. Compare hypot -> no sign issue.
// (c) threshold: solver hTh = Gamma*hypot(1, delta/Gamma)/v?? check: arnoldThresholdMT = (Gamma/v)*hypot(1,delta/Gamma) = hypot(Gamma, delta)/v == Eq.(2). OK.
// (d) 2A sideband: growth uses hx * ((P2x a) + (Q2x a*)) — BOTH P2x and Q2x enter; v uses |Q2x|/2 only.
//     Check whether P2x diagonal contributes to the threshold (it multiplies a_h, i.e. a *pump-synchronous* freq shift):
{
  for (const m of [m1, m2]) {
    const i2 = (m * M + m) * 2;
    P(`mode idx ${m}: Q2x_diag = ${Q2x[i2].toFixed(6)} + ${Q2x[i2 + 1].toFixed(6)}i (|Q| = ${(Math.hypot(Q2x[i2], Q2x[i2 + 1])).toFixed(6)}), P2x_diag = ${P2x[i2].toFixed(6)} + ${P2x[i2 + 1].toFixed(6)}i (|P| = ${(Math.hypot(P2x[i2], P2x[i2 + 1])).toFixed(6)}; ratio |P|/|Q| = ${(Math.hypot(P2x[i2], P2x[i2 + 1]) / Math.hypot(Q2x[i2], Q2x[i2 + 1])).toFixed(6)})`);
  }
}
// (e) Gamma vs doc scale: doc says Gamma=8.2-9.8 MHz at alpha=0.002. Solver at 0.0005 gives ~1/4.
P('Gamma(min)/Gamma(max) in MHz at damping=0.0005: ' + (Math.min(...gamArr) * 1000).toFixed(4) + '/' + (Math.max(...gamArr) * 1000).toFixed(4));
P('Scaled to alpha=0.002 (x4): ' + (Math.min(...gamArr) * 4000).toFixed(4) + '/' + (Math.max(...gamArr) * 4000).toFixed(4) + ' MHz — doc quotes 8.2-9.8 MHz');

fs.writeFileSync(path.join(HERE, 'nmm_check_stdout.txt'), out.join('\n') + '\n');
console.log(out.join('\n'));
