#!/usr/bin/env node
'use strict';

// Run with NODE_PATH pointing to Playwright when it is installed outside the
// checkout. SRG_CHROMIUM selects a browser executable for cloud environments.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const { chromium } = require('playwright');

const root = path.resolve(__dirname, '..');
const html = path.join(root, 'relativistic-navigator.html');
let checks = 0;
function check(condition, message) {
  assert.ok(condition, message);
  checks++;
}

async function startServer() {
  const server = http.createServer((request, response) => {
    const pathname = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
    const requested = path.resolve(root, `.${pathname === '/' ? '/relativistic-navigator.html' : pathname}`);
    if (requested !== root && !requested.startsWith(root + path.sep)) {
      response.writeHead(403).end();
      return;
    }
    fs.readFile(requested, (error, content) => {
      if (error) { response.writeHead(404).end(); return; }
      const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm', '.json': 'application/json' };
      response.setHeader('Content-Type', types[path.extname(requested)] || 'application/octet-stream');
      response.end(content);
    });
  });
  await new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(0, '127.0.0.1', resolve);
  });
  return server;
}

async function state(page, tau) {
  return page.evaluate(value => {
    const H = window.__navigatorTest;
    const current = H.seat(value);
    return { t: current.t, tau: current.tau, phi: current.phi, gamma: current.gamma, beta: current.beta, X: current.X, arrived: current.arrived };
  }, tau);
}

function validState(value, message) {
  check([value.t, value.tau, value.phi, value.gamma, value.beta, ...value.X].every(Number.isFinite), message);
  check(value.gamma >= 1 && value.beta >= 0 && value.beta <= 1, `${message}: physical speed/clock range`);
}

async function boot(context, url, content) {
  const errors = [];
  const page = await context.newPage();
  page.on('pageerror', error => errors.push(error.message));
  if (content === undefined) await page.goto(url, { waitUntil: 'load', timeout: 120000 });
  else await page.setContent(content, { waitUntil: 'load', timeout: 120000 });
  await page.waitForFunction(() => window.__navigatorTest && window.SRKitNavigator, null, { timeout: 120000 });
  await page.evaluate(() => window.__navigatorTest.stopPlay());
  check(await page.evaluate(() => !!window.SRKitNavigator.info), 'backend metadata is present');
  check(await page.evaluate(() => typeof window.SRKitNavigator.calls() === 'number' && window.SRKitNavigator.calls() > 0), 'rendering executed the Wasm backend');
  return { page, errors };
}

async function main() {
  check(fs.existsSync(html), 'build the standalone Navigator before browser testing');
  const server = await startServer();
  let browser;
  try {
    browser = await chromium.launch({
      ...(process.env.SRG_CHROMIUM ? { executablePath: process.env.SRG_CHROMIUM } : {}),   // default: Playwright's own browser
      headless: true,
      args: ['--no-sandbox', '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'],
    });
    const context = await browser.newContext({ viewport: { width: 1280, height: 900 }, reducedMotion: 'reduce' });
    await context.addInitScript(() => {
      // Keep the representative map/chart views; unrelated sky renderers are
      // tested separately from this CPU physics/backend integration smoke.
      try { localStorage.setItem('nav.panels', JSON.stringify({ chart: true, card: true, strip: true, st: false, st1: true, pano: false, look: false, plane: false, notes: false, debug: true })); } catch (_) {}
    });
    const { page, errors } = await boot(context, `http://127.0.0.1:${server.address().port}/`);
    const initialCalls = await page.evaluate(() => window.SRKitNavigator.calls());

    const destinations = await page.locator('#dests button[data-k]').evaluateAll(buttons => buttons.map(button => button.dataset.k));
    check(destinations.length > 1, 'destination controls are populated');
    for (const key of ['far', 'm31'].filter(key => destinations.includes(key))) {
      await page.locator(`#dests button[data-k="${key}"]`).click();
      check(await page.evaluate(key => window.__navigatorTest.S.sel === key, key), `destination control selects ${key}`);
      const plan = await page.evaluate(() => { const p = window.__navigatorTest.curTrip(); return { tau: p.tau, tEnd: p.tEnd }; });
      check(plan.tau > 0 && plan.tEnd >= plan.tau, `${key}: trip clocks are valid`);
      for (const tau of [-1e6, 0, plan.tau * 0.5, plan.tau, plan.tau + 1e6]) validState(await state(page, tau), `${key}: state at ship time ${tau}`);
    }
    await page.locator('#tlaunch').click();
    check(await page.evaluate(() => window.__navigatorTest.S.tau === 0), 'Launch returns to the start');
    await page.locator('#tstepback button').last().click();
    check((await page.locator('#tstepnote').innerText()).includes('before launch'), 'backwards step renders the waiting state');
    await page.locator('#tarrive').click();
    await page.locator('#tstepfwd button').first().click();
    check((await page.locator('#tstepnote').innerText()).includes('after'), 'forward step renders the post-arrival state');
    await page.locator('#tflip').click();
    validState(await page.evaluate(() => { const H = window.__navigatorTest, s = H.seat(H.S.tau); return { ...s, X: s.X }; }), 'Flip updates the actual ride');

    // Exercise the four chart controls without relying on timing budgets or
    // screenshot equality; all receive their geometry from the selected ride.
    await page.locator('#s1stars').selectOption('0');
    await page.locator('#s1grid').selectOption('off');
    for (const chart of ['gal', 'mcrf', 'lc', 'radar']) {
      const control = page.locator(`#s1chart button[data-c="${chart}"]`);
      await control.click();
      check(await control.getAttribute('aria-pressed') === 'true', `${chart} chart control is active`);
    }

    await page.locator('#mode button[data-m="loop"]').click();
    check(await page.evaluate(() => window.__navigatorTest.S.mode === 'loop'), 'circular-orbit control switches mode');
    const orbit = await page.evaluate(() => { const L = window.__navigatorTest.loopP(); return { lapTau: L.lapTau, gamma: L.gamma, beta: L.beta }; });
    check(orbit.lapTau > 0 && orbit.gamma > 1 && orbit.beta > 0 && orbit.beta < 1, 'orbit parameters are physical');
    for (const tau of [-orbit.lapTau, 0, orbit.lapTau * 0.25, orbit.lapTau]) validState(await state(page, tau), `orbit state at ${tau}`);
    await page.locator('#loopr').selectOption('1');
    check(await page.evaluate(() => window.__navigatorTest.S.loopR === 1), 'ring selector changes the orbit');

    check(!/\bNaN\b/.test(await page.locator('#facts').innerText()), 'trip/orbit card contains no NaN');
    check(!/\bNaN\b/.test(await page.locator('#readout').innerText()), 'clock readouts contain no NaN');
    check(await page.evaluate(before => window.SRKitNavigator.calls() > before, initialCalls), 'interactive changes continue executing Wasm');
    check(errors.length === 0, `no desktop page errors: ${errors.join('; ')}`);

    await page.setViewportSize({ width: 390, height: 844 });
    await page.waitForTimeout(250);
    check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), '390 px layout has no horizontal overflow');
    await page.locator('#mode button[data-m="transit"]').click();
    await page.locator('#tlaunch').click();
    check(await page.locator('#tstepnote').isVisible(), 'phone ship-time controls remain visible');
    check(errors.length === 0, `no phone page errors: ${errors.join('; ')}`);
    await context.close();

    const offline = await browser.newContext({ viewport: { width: 1000, height: 800 }, reducedMotion: 'reduce' });
    const network = [];
    await offline.route(/^https?:\/\//, route => { network.push(route.request().url()); return route.abort(); });
    let file;
    let directFileNavigation = true;
    try {
      file = await boot(offline, pathToFileURL(html).href);
    } catch (error) {
      if (!String(error).includes('ERR_BLOCKED_BY_ADMINISTRATOR')) throw error;
      directFileNavigation = false;
      console.log('SKIP direct file:// navigation: managed Chromium policy blocks the file scheme; checking the same standalone HTML offline with setContent.');
      file = await boot(offline, null, fs.readFileSync(html, 'utf8'));
    }
    validState(await state(file.page, 1), 'standalone file opens with valid Wasm state');
    check(file.errors.length === 0, `no standalone file page errors: ${file.errors.join('; ')}`);
    check(!network.some(url => /\.wasm(?:\?|$)/.test(url)), 'standalone backend needs no Wasm network fetch');
    await offline.close();
    console.log(`Navigator browser smoke: ${checks} checks passed (HTTP desktop, phone and offline standalone HTML; direct file navigation ${directFileNavigation ? 'passed' : 'blocked by browser policy'}).`);
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
}

main().catch(error => { console.error(error.stack || error); process.exitCode = 1; });
