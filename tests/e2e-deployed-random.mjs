import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';

const args = process.argv.slice(2);
const argVal = (flag, fb) => {
  const hit = args.find(a => a.startsWith(flag + '='));
  return hit ? hit.slice(flag.length + 1) : fb;
};

const TARGET_URL = argVal('--url', process.env.TARGET_URL || 'https://kfwbnh5y25-source.github.io/spacetime-grid-/spin-waves/');
const SEED = Number(argVal('--seed', process.env.E2E_SEED || ((Date.now() & 0x7fffffff) >>> 0)));
const SHOT_DIR = argVal('--shots', '');

function makeRng(seed) {
  let s = (seed >>> 0) || 1;
  return () => {
    let t = (s = (s + 0x6D2B79F5) >>> 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const rand = makeRng(SEED);
const randBetween = (lo, hi) => lo + rand() * (hi - lo);
const pick = arr => arr[Math.floor(rand() * arr.length)];

function b64ToF32(b64) {
  const buf = Buffer.from(b64, 'base64');
  const u8 = new Uint8Array(buf.length);
  u8.set(buf);
  return new Float32Array(u8.buffer);
}

// Extract ground-truth NMM_DATA directly from served HTML (independent analytical oracle)
async function fetchGroundTruthOracle(url) {
  const resp = await fetch(url, { cache: 'no-store' });
  assert.ok(resp.ok, `GET ${url} failed with HTTP ${resp.status}`);
  const html = await resp.text();
  const m = html.match(/const NMM_DATA=(\{[\s\S]*?\});/);
  assert.ok(m, 'Could not locate NMM_DATA JSON in served HTML');
  const d = JSON.parse(m[1]);
  const M = d.M;
  const F_M = d.fM;
  const MU0_MS_MT = d.mu0Ms_mT;
  const fGhz = b64ToF32(d.fGhz);
  const Kp = b64ToF32(d.Kp);
  const Q2x = b64ToF32(d.Q2x);
  const N4A = b64ToF32(d.N4A);

  const gammaGHz = (mIdx, alpha) => Math.max(1e-7, alpha * Kp[(mIdx * M + mIdx) * 2] * F_M);
  const vParamGHz_mT = mIdx => 0.5 * Math.hypot(Q2x[(mIdx * M + mIdx) * 2], Q2x[(mIdx * M + mIdx) * 2 + 1]) * F_M / MU0_MS_MT;
  const hThMT = (mIdx, alpha) => gammaGHz(mIdx, alpha) / Math.max(1e-9, vParamGHz_mT(mIdx));
  const selfNfsGHz = mIdx => N4A[mIdx * M + mIdx] * F_M;
  const steadyStateAmpEq33 = (mIdx, fPumpGHz, hPumpMT, alpha) => {
    const gam = gammaGHz(mIdx, alpha);
    const vm = vParamGHz_mT(mIdx);
    if (hPumpMT * vm <= gam) return 0;
    const rad = Math.sqrt((vm * hPumpMT) ** 2 - gam * gam);
    const delta = 0.5 * (fPumpGHz - 2 * fGhz[mIdx]);
    if (Math.abs(delta) >= rad) return 0;
    const nkk = selfNfsGHz(mIdx);
    const sgn = nkk >= 0 ? 1 : -1;
    if (sgn * delta < -rad) return 0;
    const a2 = (delta + sgn * rad) / nkk;
    return a2 > 0 ? Math.sqrt(a2) : 0;
  };
  return { M, F_M, MU0_MS_MT, fGhz, selfNfsGHz, gammaGHz, vParamGHz_mT, hThMT, steadyStateAmpEq33 };
}

const chromeBin = ['/usr/bin/chromium', '/usr/bin/google-chrome', '/usr/bin/chromium-browser'].find(p => fs.existsSync(p));
assert.ok(chromeBin, 'Chromium binary not found on VM');

const userDataDir = `/tmp/e2e-prof-${process.pid}-${SEED}`;
fs.rmSync(userDataDir, { recursive: true, force: true });

console.log(`[E2E] Target: ${TARGET_URL} (seed=${SEED})`);
const oracle = await fetchGroundTruthOracle(TARGET_URL);

const chrome = spawn(chromeBin, [
  '--headless=new', '--no-sandbox', '--disable-gpu', '--remote-debugging-port=0',
  `--user-data-dir=${userDataDir}`, '--window-size=1440,900',
  TARGET_URL
], { stdio: ['ignore', 'pipe', 'pipe'] });

let ws = null;
try {
  const wsUrl = await new Promise((resolve, reject) => {
    let buf = '';
    const onData = d => {
      buf += d.toString();
      const match = buf.match(/DevTools listening on (ws:\/\/[^\s]+)/);
      if (match) resolve(match[1]);
    };
    chrome.stderr.on('data', onData);
    chrome.stdout.on('data', onData);
    setTimeout(() => reject(new Error('Chromium CDP timeout: ' + buf)), 15000);
  });

  const targets = await (await fetch(wsUrl.replace(/^ws:/, 'http:').replace(/\/devtools\/browser\/.*/, '/json/list'))).json();
  const pageTarget = targets.find(t => t.type === 'page') || targets[0];
  assert.ok(pageTarget, 'No Chromium page target opened');

  ws = new WebSocket(pageTarget.webSocketDebuggerUrl);
  await new Promise(r => ws.addEventListener('open', r, { once: true }));

  let msgId = 0;
  const pending = new Map();
  const runtimeErrors = [];
  const consoleErrors = [];
  ws.addEventListener('message', ev => {
    const msg = JSON.parse(ev.data);
    if (msg.id && pending.has(msg.id)) {
      const { resolve, reject } = pending.get(msg.id);
      pending.delete(msg.id);
      if (msg.error) reject(new Error(JSON.stringify(msg.error)));
      else resolve(msg.result);
    } else if (msg.method === 'Runtime.exceptionThrown') {
      runtimeErrors.push(JSON.stringify(msg.params.exceptionDetails));
    } else if (msg.method === 'Runtime.consoleAPICalled' && msg.params.type === 'error') {
      const txt = msg.params.args.map(a => a.value || a.description || '').join(' ');
      if (!txt.includes('401') && !txt.includes('Failed to load resource')) consoleErrors.push(txt);
    }
  });

  const cdp = (method, params = {}) => new Promise((resolve, reject) => {
    const id = ++msgId;
    pending.set(id, { resolve, reject });
    ws.send(JSON.stringify({ id, method, params }));
  });
  const evalJs = async expr => {
    const r = await cdp('Runtime.evaluate', { expression: expr, returnByValue: true, awaitPromise: true });
    if (r?.exceptionDetails) throw new Error(`Page eval error: ${JSON.stringify(r.exceptionDetails)}`);
    return r?.result?.value;
  };
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const clickId = async id => {
    const ok = await evalJs(`(() => { const el = document.getElementById(${JSON.stringify(id)}); if (!el) return false; el.click(); return true; })()`);
    assert.equal(ok, true, `Element #${id} not found for click`);
  };
  const setInputVal = async (id, val, evType = 'input') => {
    const ok = await evalJs(`(() => {
      const el = document.getElementById(${JSON.stringify(id)});
      if (!el) return false;
      el.value = ${JSON.stringify(String(val))};
      el.dispatchEvent(new Event(${JSON.stringify(evType)}, { bubbles: true }));
      return true;
    })()`);
    assert.equal(ok, true, `Input #${id} not found`);
  };

  await cdp('Page.enable');
  await cdp('Runtime.enable');
  await cdp('Emulation.setDeviceMetricsOverride', { width: 1440, height: 900, deviceScaleFactor: 1, mobile: false });

  // Wait for page bootstrap
  let loaded = false;
  for (let i = 0; i < 60; i++) {
    const txt = await evalJs(`document.getElementById('domModeBadge')?.textContent || ''`);
    if (txt.includes('MODE:')) { loaded = true; break; }
    await sleep(250);
  }
  assert.ok(loaded, 'Page did not finish initializing #domModeBadge');

  // 1. Fresh What's New modal + Tutorial close + Release check via real UI clicks
  const wnVisible = await evalJs(`!document.getElementById('whatsNewModal').hidden`);
  assert.equal(wnVisible, true, 'Fresh launch must show #whatsNewModal');
  await clickId('wnGotItBtn');
  assert.equal(await evalJs(`document.getElementById('whatsNewModal').hidden`), true);
  const seenVer = await evalJs(`localStorage.getItem('spin-waves-whatsnew-seen-version')`);
  assert.equal(seenVer, 'spacetime-shell-v84', `Expected deployed version spacetime-shell-v84, got ${seenVer}`);
  await clickId('swTutClose');
  assert.equal(await evalJs(`document.getElementById('swTutorial').hidden`), true);

  // Helper to run FAST->SS via real UI button and wait until settled
  const runFastForwardViaUi = async (timeoutMs = 45000) => {
    await evalJs(`(() => {
      const q = document.getElementById('quickFastBtn');
      if (q && !q.textContent.includes('STOP')) q.click();
    })()`);
    const tStart = Date.now();
    while (Date.now() - tStart < timeoutMs) {
      await sleep(200);
      const txt = await evalJs(`document.getElementById('quickFastBtn')?.textContent || ''`);
      if (txt === 'FAST→SS') return;
    }
    throw new Error('FAST->SS did not settle within timeout');
  };

  // Inspect Live Calibration via real modal button clicks and return #vmLiveOut text
  const inspectLiveViaUi = async () => {
    await clickId('quickVerifyBtn');
    await evalJs(`document.querySelector('#vmTabs .vmTab[data-tab="live"]')?.click()`);
    await clickId('vmInspectLiveBtn');
    const out = await evalJs(`document.getElementById('vmLiveOut')?.textContent || ''`);
    await clickId('vmCloseBtn');
    return out;
  };

  // 2. Randomized Single-Mode Eq. (33) UI Trials (Tab 1)
  await clickId('btnNmmRwa');
  await clickId('driveParallel');
  await setInputVal('slTemp', 0, 'input');

  const trialModes = [13, 14, 16, 19, 20, 27, 28];
  for (let t = 0; t < 3; t++) {
    const m = pick(trialModes);
    const alphaReq = +randBetween(0.0014, 0.0028).toFixed(4);
    const ratio = randBetween(1.22, 1.48);
    const hThReq = oracle.hThMT(m, alphaReq);
    const hPumpReq = Math.min(2.45, +(ratio * hThReq).toFixed(2));
    // Place detuning strictly on the active side of the tongue fold to avoid bistability ambiguity
    const nkkSign = oracle.selfNfsGHz(m) >= 0 ? 1 : -1;
    const fPumpReq = +(2 * oracle.fGhz[m] + nkkSign * 0.0015).toFixed(4);

    await setInputVal('slDamp', alphaReq, 'input');
    await setInputVal('slFreq', fPumpReq, 'input');
    await setInputVal('slAmp', hPumpReq, 'input');
    const alpha = await evalJs(`parseFloat(document.getElementById('slDamp').value)`);
    const fPump = await evalJs(`parseFloat(document.getElementById('slFreq').value)`);
    const hPump = await evalJs(`parseFloat(document.getElementById('slAmp').value)`);
    const hTh = oracle.hThMT(m, alpha);
    const expAmp = oracle.steadyStateAmpEq33(m, fPump, hPump, alpha);
    const expShiftMHz = oracle.selfNfsGHz(m) * expAmp * expAmp * 1000;

    await clickId('btnResetTime');
    await runFastForwardViaUi();

    const liveTxt = await inspectLiveViaUi();
    const threshTxt = await evalJs(`document.getElementById('threshBadge')?.textContent || ''`);
    const reMode = new RegExp(`m${m}\\([^)]*\\ba=([0-9.]+),\\s*h_th=([0-9.]+)mT,\\s*Δ4A=([-0-9.]+)MHz`);
    const mm = liveTxt.match(reMode);
    assert.ok(mm, `Trial ${t + 1}: m${m} not found active in #vmLiveOut:\n${liveTxt}`);
    const simAmp = parseFloat(mm[1]);
    const simHth = parseFloat(mm[2]);
    const simShift = parseFloat(mm[3]);
    assert.ok(Math.abs(simHth - hTh) < 0.005, `Trial ${t + 1} m${m} h_th mismatch: UI=${simHth} vs oracle=${hTh.toFixed(4)}`);
    assert.ok(Math.abs(simAmp - expAmp) < Math.max(0.012, expAmp * 0.08),
      `Trial ${t + 1} m${m} (alpha=${alpha}, f=${fPump}G, h=${hPump}mT) steady-state amp mismatch: UI=${simAmp} vs Eq33=${expAmp.toFixed(5)}`);
    assert.ok(Math.sign(simShift) === Math.sign(expShiftMHz) || Math.abs(expShiftMHz) < 0.5,
      `Trial ${t + 1} m${m} 4A sign mismatch: UI=${simShift} vs oracle=${expShiftMHz.toFixed(2)}`);
    assert.ok(threshTxt.includes('ABOVE'), `Trial ${t + 1} m${m} expected ABOVE threshold badge, got ${threshTxt}`);
    const arnoldNote = await evalJs(`document.getElementById('arnoldTermsNote')?.textContent || ''`);
    assert.ok(arnoldNote.includes('4A self:'), `Trial ${t + 1} m${m} expected single-mode '4A self:' readout on Arnold plot, got: ${arnoldNote}`);
    if (SHOT_DIR && t === 0) {
      fs.mkdirSync(SHOT_DIR, { recursive: true });
      const sSm = await cdp('Page.captureScreenshot', { format: 'png' });
      fs.writeFileSync(path.join(SHOT_DIR, 'e2e_tab1_singlemode.png'), Buffer.from(sSm.data, 'base64'));
    }
    console.log(`[E2E] Single-Mode Trial ${t + 1}: m${m} α=${alpha} h_p=${hPump}mT -> a_sim=${simAmp.toFixed(5)} vs Eq33=${expAmp.toFixed(5)} (${arnoldNote.split('·')[0].trim()}) PASS`);
  }

  // 3. Randomized 50-Ohm RF Calculator UI Trials
  for (let t = 0; t < 3; t++) {
    const vpk = Math.round(randBetween(120, 780));
    const gainDB = +randBetween(-10, 10).toFixed(1);
    await setInputVal('rfVpkTargetIn', vpk, 'input');
    await setInputVal('rfChainGainIn', gainDB, 'input');
    const outTxt = await evalJs(`document.getElementById('rfTargetOut')?.textContent || ''`);
    const expMW = (vpk * vpk) / 100000;
    const expDBm = 10 * Math.log10(expMW);
    const m = outTxt.match(/=\s*([0-9.]+)\s*mW\s*\(([-0-9.]+)\s*dBm\)/);
    assert.ok(m, `RF output parse failed: ${outTxt}`);
    assert.ok(Math.abs(parseFloat(m[1]) - expMW) < 0.005, `P_AWG mW mismatch: ${m[1]} vs ${expMW.toFixed(3)}`);
    assert.ok(Math.abs(parseFloat(m[2]) - expDBm) < 0.03, `P_AWG dBm mismatch: ${m[2]} vs ${expDBm.toFixed(2)}`);
  }
  console.log('[E2E] RF 50-Ohm randomized UI conversions PASS');

  // 4. Dual-Tone Scenario 1 Coex & Scenario 3 4B Phase Gate via UI Clicks
  await clickId('quickCoexBtn');
  await runFastForwardViaUi();
  const coexTxt = await inspectLiveViaUi();
  assert.ok(coexTxt.includes('m13(') && coexTxt.includes('m23(') && coexTxt.includes('"scenario":1'),
    `Coex scenario failed:\n${coexTxt}`);

  await clickId('presetCross4B');
  await setInputVal('slPhase1', 0, 'input');
  await runFastForwardViaUi();
  const c4b0 = await inspectLiveViaUi();
  const m15_0 = c4b0.match(/m15\([^)]*\ba=([0-9.]+)/);
  assert.ok(m15_0 && parseFloat(m15_0[1]) > 0.010, `4B φ1=0° expected m15 ignition, got:\n${c4b0}`);

  await clickId('presetCross4B');
  await setInputVal('slPhase1', 180, 'input');
  await runFastForwardViaUi();
  const c4b180 = await inspectLiveViaUi();
  const m15_180 = c4b180.match(/m15\([^)]*\ba=([0-9.]+)/);
  assert.ok(!m15_180 || parseFloat(m15_180[1]) < 0.003, `4B φ1=180° expected m15 suppression, got:\n${c4b180}`);
  console.log('[E2E] Dual-Tone Coex & 4B Phase Gate UI checks PASS');

  if (SHOT_DIR) {
    fs.mkdirSync(SHOT_DIR, { recursive: true });
    await clickId('quickCoexBtn');
    await sleep(600);
    const s1 = await cdp('Page.captureScreenshot', { format: 'png' });
    fs.writeFileSync(path.join(SHOT_DIR, 'e2e_tab1_desktop.png'), Buffer.from(s1.data, 'base64'));
  }

  // 5. Tab 2 (2 · ARNOLD SWEEP) Randomized UI + Canvas Click + Send-to-Tab-1 Verification
  await clickId('tabSweepBtn');
  assert.equal(await evalJs(`document.getElementById('dashShell').hidden`), true);
  assert.equal(await evalJs(`document.getElementById('sweepShell').hidden`), false);

  await clickId('sweepPresetZoom');
  assert.equal(await evalJs(`Number(document.getElementById('sweepFMinIn').value)`), 4050);
  assert.equal(await evalJs(`Number(document.getElementById('sweepFMaxIn').value)`), 4080);

  await clickId('sweepPresetWide');
  const rKEff = +randBetween(1.30, 1.80).toFixed(2);
  const rDamp = +randBetween(0.0016, 0.0026).toFixed(4);
  const rTemp = rand() < 0.5 ? 0 : 300;
  await setInputVal('sweepKEffIn', rKEff, 'change');
  await setInputVal('sweepDampIn', rDamp, 'change');
  await setInputVal('sweepTempIn', rTemp, 'change');
  await clickId('sweepRunFastBtn');

  await clickId('sweepMetricWidth');
  await clickId('sweepMetricNeff');
  await clickId('sweepMetricSignal');

  // Click directly on the Kittel m27 tongue coordinate (4109.2 MHz, 450 mV_pk) on #sweepCanvas
  const clickPt = await evalJs(`(() => {
    const c = document.getElementById('sweepCanvas');
    const r = c.getBoundingClientRect();
    const padL = 50, padR = 14, padT = 16, padB = 34;
    const pw = Math.max(40, r.width - padL - padR), ph = Math.max(40, r.height - padT - padB);
    const fx = (4109.2 - 3760) / (4250 - 3760);
    const fy = (450 - 200) / (500 - 200);
    const cx = r.left + padL + fx * pw;
    const cy = r.top + padT + (1 - fy) * ph;
    c.dispatchEvent(new PointerEvent('pointerdown', { clientX: cx, clientY: cy, bubbles: true }));
    return { hover: document.getElementById('sweepHoverInfo').textContent, slice: document.getElementById('sweepSliceMeta').textContent };
  })()`);
  const hpM = clickPt.hover.match(/μ₀h_p=([0-9.]+)\s*mT/);
  const sigM = clickPt.hover.match(/Sig=([-0-9.]+)\s*dB\s*\(SNR=([0-9.]+)×\)/);
  assert.ok(hpM && sigM, `Sweep hover info parse failed: ${clickPt.hover}`);
  const expHp = (450 * 1e-3) * rKEff;
  assert.ok(Math.abs(parseFloat(hpM[1]) - expHp) < 0.03, `Sweep h_p mismatch: ${hpM[1]} vs ${expHp.toFixed(3)}`);
  assert.ok(parseFloat(sigM[2]) >= 10, `Kittel m27 point (4109.2 MHz, 450 mV, k_eff=${rKEff}) expected SNR>=10x, got ${sigM[2]}`);

  if (SHOT_DIR) {
    await sleep(300);
    const s2 = await cdp('Page.captureScreenshot', { format: 'png' });
    fs.writeFileSync(path.join(SHOT_DIR, 'e2e_tab2_desktop.png'), Buffer.from(s2.data, 'base64'));
  }

  // Send selected Sweep point back into Tab 1 and verify exact transfer
  await clickId('sweepSendLiveBtn');
  assert.equal(await evalJs(`document.getElementById('dashShell').hidden`), false);
  assert.equal(await evalJs(`document.getElementById('sweepShell').hidden`), true);
  const transferredFreq = await evalJs(`parseFloat(document.getElementById('slFreq').value)`);
  const transferredAmp = await evalJs(`parseFloat(document.getElementById('slAmp').value)`);
  assert.ok(Math.abs(transferredFreq - 4.109) < 0.015, `Sent freq mismatch: ${transferredFreq}`);
  assert.ok(Math.abs(transferredAmp - expHp) < 0.04, `Sent h_p mismatch: ${transferredAmp} vs ${expHp.toFixed(3)}`);
  console.log(`[E2E] Tab 2 Sweep randomized k_eff=${rKEff}, α=${rDamp}, T=${rTemp}K & Send->Tab 1 PASS`);

  // 6. Run full 8-Stage Verification Protocol via real UI click
  await clickId('quickVerifyBtn');
  await clickId('vmRunAllBtn');
  let overallBadge = '', vmStatus = '';
  for (let i = 0; i < 180; i++) {
    await sleep(500);
    vmStatus = await evalJs(`document.getElementById('vmStatusText')?.textContent || ''`);
    overallBadge = await evalJs(`document.getElementById('vmScoreBadge')?.textContent || ''`);
    if (vmStatus.includes('Completed 8/8')) break;
  }
  assert.ok(vmStatus.includes('Completed 8/8'), `8-stage verify did not complete 8/8 stages: ${vmStatus} (${overallBadge})`);
  assert.ok(overallBadge.includes('PASS') && !overallBadge.includes('FAIL'), `Deployed 8-stage verify failed: ${overallBadge}`);
  await clickId('vmCloseBtn');
  console.log(`[E2E] Deployed 8-Stage Verification Protocol UI click -> ${overallBadge} PASS`);

  // 7. Mobile (390x844) responsiveness check + screenshot
  await cdp('Emulation.setDeviceMetricsOverride', { width: 390, height: 844, deviceScaleFactor: 2, mobile: true });
  await clickId('tabSweepBtn');
  await sleep(350);
  const mobBox = await evalJs(`(() => {
    const r = document.getElementById('sweepCanvas').getBoundingClientRect();
    const mapR = document.getElementById('sweepMapCard').getBoundingClientRect();
    const sliceR = document.getElementById('sweepSliceCard').getBoundingClientRect();
    return { w: r.width, h: r.height, scrollW: document.documentElement.scrollWidth, mapBottom: mapR.bottom, sliceTop: sliceR.top };
  })()`);
  assert.ok(mobBox.w >= 260 && mobBox.h >= 200, `Mobile Tab 2 canvas degenerate: ${JSON.stringify(mobBox)}`);
  assert.ok(mobBox.mapBottom <= mobBox.sliceTop + 1, `Mobile Tab 2 card overlap detected: mapBottom=${mobBox.mapBottom} > sliceTop=${mobBox.sliceTop}`);
  assert.ok(mobBox.scrollW <= 410, `Mobile horizontal overflow detected: scrollWidth=${mobBox.scrollW}`);
  if (SHOT_DIR) {
    const s3 = await cdp('Page.captureScreenshot', { format: 'png' });
    fs.writeFileSync(path.join(SHOT_DIR, 'e2e_tab2_mobile.png'), Buffer.from(s3.data, 'base64'));
  }

  assert.deepEqual(runtimeErrors, [], `Runtime page exceptions: ${runtimeErrors.join('; ')}`);
  assert.deepEqual(consoleErrors, [], `Console errors: ${consoleErrors.join('; ')}`);
  console.log('[E2E] ALL DEPLOYED RANDOMIZED UI + ANALYTICAL VERIFICATION CHECKS PASSED!');
} finally {
  try { ws?.close(); } catch (_) {}
  try { chrome?.kill('SIGKILL'); } catch (_) {}
  await new Promise(r => setTimeout(r, 150));
  try { fs.rmSync(userDataDir, { recursive: true, force: true }); } catch (_) {}
}
