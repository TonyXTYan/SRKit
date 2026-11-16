// Tests for SRKit's JavaScript binding (dist/srkit.js), run in Node:
//     python3 -I tools/build_js.py && node --expose-gc tools/test_js.js
// They check the binding, not the physics (the Rust tests do that): that every part of the API reaches the module
// and comes back in the documented shape, that plain-object worldlines give the same answers whether the module
// holds a copy or not, that memory growth and garbage collection are handled, and that failures are explicit.
// Exits non-zero on any failure. Code AI-generated (Claude).
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs'), path = require('node:path'), vm = require('node:vm');
const dist = path.resolve(__dirname, '../dist');
const SRKit = require(path.join(dist, 'srkit.js'));
let checks = 0;
const ok = (cond, what) => { assert.ok(cond, what); checks++; };
const near = (a, b, what, rel = 1e-12) => { ok(Number.isFinite(a) && Math.abs(a - b) <= rel * Math.max(1, Math.abs(a), Math.abs(b)), `${what}: ${a} against ${b}`); };
const g = SRKit.units.G_PER_YEAR, S = SRKit.seg, mk = S.mk;

// ---- the module and its metadata ----
ok(SRKit.info.backend === 'srkit-wasm' && /^\d+\.\d+\.\d+$/.test(SRKit.info.version) && /^[0-9a-f]{64}$/.test(SRKit.info.sha256), 'info names the build');
const wasm = fs.readFileSync(path.join(dist, 'srkit.wasm'));
ok(require('node:crypto').createHash('sha256').update(wasm).digest('hex') === SRKit.info.sha256 && wasm.length === SRKit.info.bytes, 'the embedded module is dist/srkit.wasm, by SHA-256');
near(SRKit.units.G_PER_YEAR, 9.80665 * 31557600 / 299792458, '1 g per year', 0);
near(SRKit.units.LY, 9460730472580800, 'the light-year is exact', 0);

// ---- segments of every kind ----
const kinds = {
  rest: mk.rest(0, 10, [1, 2, 3]),
  lin: mk.lin(0, 10, [1, 2, 3], [0.6, 0, 0]),
  hyp: mk.hyp(0, 10, [0, 0, 0], [0, 1, 0], 1, 0),
  circ: mk.circ(0, 10, [0, 0, 0], 2, 0.25, [1, 0, 0], [0, 1, 0], 0),
  null: mk.null(0, 5, [0, 0, 0], [0, 0, 1]),
  steer: mk.steer(0, 10, [0, 0, 0], [0.75, 0, 0], [0, 1, 0], 1),
  burn: mk.burn(0, 0.5, [0, 0, 0], [1, 0, 0], 0, 1, 1),
};
for (const [kind, seg] of Object.entries(kinds)) {
  const st = S.at(seg, 0.4), s = S.dtau(seg, 0.4);
  ok(st.X.length === 3 && st.V.length === 3 && st.X.every(Number.isFinite), `${kind}: at() gives a position and a velocity`);
  ok(kind === 'null' ? s === 0 : s > 0 && s <= 0.4, `${kind}: proper time is no more than coordinate time`);
  if (kind !== 'null') { near(S.tOf(seg, s), 0.4, `${kind}: tOf inverts dtau`, 1e-12); const e = S.atTau(seg, s); near(e.t, 0.4, `${kind}: atTau agrees`, 1e-12); ok(Math.hypot(...e.X.map((c, i) => c - st.X[i])) < 1e-12, `${kind}: atTau position`); }
  ok(S.vmax(seg) >= 0 && S.vmax(seg) <= 1, `${kind}: vmax is a fraction of c`);
  const tr = S.recv(seg, -3, [5, 5, 0], -Infinity, Infinity);
  if (tr != null) { const p = S.at(seg, tr).X; near(tr - -3, Math.hypot(p[0] - 5, p[1] - 5, p[2]), `${kind}: recv is on the light cone`, 1e-11); }
  const ts = S.send(seg, 30, [5, 5, 0], -Infinity, Infinity);
  if (ts != null) { const p = S.at(seg, ts).X; near(30 - ts, Math.hypot(p[0] - 5, p[1] - 5, p[2]), `${kind}: send is on the light cone`, 1e-11); }
}
near(S.at(kinds.hyp, 1).X[1], Math.SQRT2 - 1, 'hyperbolic motion from rest: x = sqrt(1 + t^2) - 1');
near(S.dtau(kinds.hyp, 1), Math.asinh(1), 'hyperbolic motion: tau = asinh(t)');
ok(S.atTau(kinds.null, 0).jump === true && S.atTau(kinds.null, 0).t === 5, 'a jump takes no proper time and lands at its end');
// thrust across the motion: the speed along the old heading falls as gamma grows, with no force along it
ok(S.at(kinds.steer, 3).V[0] < 0.6 && S.at(kinds.steer, 3).V[1] > 0.5, 'steered thrust turns the ship');
near(S.atTau(kinds.burn, 0.5).phi, -Math.log(1 - 0.5), 'a photon rocket that has burned half its mass has rapidity ln 2');
// a horizon: light sent from the launch point after t = 1/a never catches a ship that keeps accelerating
const fleeing = mk.hyp(0, Infinity, [0, 0, 0], [1, 0, 0], 2, 0);
ok(S.recv(fleeing, 0.4, [0, 0, 0], 0, Infinity) > 0.4 && S.recv(fleeing, 0.6, [0, 0, 0], 0, Infinity) === null, 'recv returns null behind the horizon');
ok(S.horizonMargin(fleeing, 0.4, [0, 0, 0]) > 0 && S.horizonMargin(fleeing, 0.6, [0, 0, 0]) < 0 && S.horizonMargin(kinds.lin, 0, [0, 0, 0]) === null, 'horizonMargin');
assert.throws(() => S.at({ kind: 'warp', t0: 0, t1: 1 }, 0), /segment kind/); checks++;

// ---- a plan, and its worldline three ways ----
const start = { t: 2100, tau: 0, X: [0, 0, 0] }, proxima = [4.2465, 0, 0];
const ride = SRKit.plan.line(start, proxima, g);
ok(ride.pieces.length === 2 && ride.pieces[0].phase === 'boost' && ride.pieces[1].phase === 'brake', 'a rest-to-rest run boosts then brakes');
near(ride.end.tau, 3.5421, 'years aboard to Proxima at 1 g', 1e-4); near(ride.end.t - 2100, 5.8726, 'years on Earth', 1e-4);
ok(ride.end.X[0] === 4.2465 && ride.dphi > 3.6 && ride.coastL === 0, 'the run ends on its target');
const capped = SRKit.plan.line(start, proxima, g, 0.5);
ok(capped.pieces.map(p => p.phase).join() === 'boost,coast,brake' && capped.coastL > 0 && capped.end.tau > ride.end.tau, 'a capped run coasts');
ok(SRKit.plan.line(start, [1e-3, 0, 0], g, Infinity, 0, 2) === null, 'a run too short to reach the speed asked at its end is refused with null');
const segs = () => [mk.rest(2000, 2100, [0, 0, 0]), ...SRKit.plan.line(start, proxima, g).pieces.map(p => p.g), mk.rest(ride.end.t, Infinity, proxima)];
const W = SRKit.wl.make(segs(), 2100, 0), hand = { segs: W.segs }, owned = SRKit.Worldline.from(segs(), 2100, 0);
ok(W.tau0s.length === 4 && W.tau0s[0] === -100 && W.tau0s[1] === 0 && Math.abs(W.tau0s[3] - ride.end.tau) < 1e-12 && W.hasNull === false, 'make() anchors proper time');
ok(owned.continuous() && owned.length === 4, 'the worldline joins');
SRKit.wl.options.handles = false; const plain = SRKit.wl.make(segs(), 2100, 0); SRKit.wl.options.handles = true;
for (let k = 0; k < 200; k++) {
  const t = 2095 + 0.07 * k, X = [3 * Math.sin(k), 2 * Math.cos(k), 0.5];
  const a = SRKit.wl.at(W, t), b = SRKit.wl.at(hand, t), c = owned.at(t), d = SRKit.wl.at(plain, t);
  ok(a.X.every((v, i) => v === b.X[i] && v === c.X[i] && v === d.X[i]) && a.V.every((v, i) => v === b.V[i] && v === c.V[i]), 'at(): the copy the module holds, the object and the owned worldline agree to the bit');
  ok(SRKit.wl.tau(W, t) === SRKit.wl.tau(hand, t) && SRKit.wl.tau(W, t) === c.tau && owned.tau(t) === c.tau, 'tau() agrees');
  ok(SRKit.wl.recv(W, t, X) === SRKit.wl.recv(hand, t, X) && SRKit.wl.recv(W, t, X) === owned.recv(t, X), 'recv() agrees');
  ok(SRKit.wl.send(W, t + 9, X) === SRKit.wl.send(hand, t + 9, X) && SRKit.wl.send(W, t + 9, X) === owned.send(t + 9, X), 'send() agrees');
}
const seen = SRKit.wl.send(W, 2105, [0, 0, 0]);
near(2105 - seen, SRKit.wl.at(W, seen).X[0], 'what Earth sees in 2105 left the ship one light-distance earlier', 1e-12);
near(owned.tOfTau(1.5), SRKit.seg.tOf(W.segs[1], 1.5), 'tOfTau');
ok(owned.find(2101) === 1 && SRKit.wl.find(W, 2101) === 1 && SRKit.wl.findTau(W, 2) === 2, 'find');
const radar = owned.radar(2104, [1, 0, 0]); ok(radar.tauSend < radar.tauRecv, 'radar brackets the event');
const back = owned.segment(1); ok(back.kind === 'hyp' && back.al === g && back.tau0 === 0 && owned.segment(9) === null, 'segment() reads a segment back');
// extended after make(): the stale copy is not used
const grown = SRKit.wl.make(segs(), 2100, 0); grown.segs.push(mk.rest(3000, Infinity, [9, 9, 9])); grown.t0s.push(3000);
ok(SRKit.wl.at(grown, 3001).X[0] === 9, 'a worldline changed after make() is answered from the object, not the stale copy');
// batches
const ts = Float64Array.from({ length: 1000 }, (_, i) => 2099 + i * 0.01), rows = owned.atBatch(ts);
ok(rows.length === 7000 && [0, 137, 999].every(i => { const c = owned.at(ts[i]); return rows[7 * i] === c.X[0] && rows[7 * i + 3] === c.V[0] && rows[7 * i + 6] === c.tau; }), 'atBatch equals at(), row for row');
const events = []; for (let i = 0; i < 300; i++) events.push(2110 + i * 0.01, 0, 0, 0);
const cone = owned.coneBatch(events, true);
ok(cone.length === 300 && cone[17] === owned.send(events[68], [0, 0, 0]), 'coneBatch equals send()');
owned.free(); ok(Number.isNaN(owned.tau(2100)), 'a freed worldline answers NaN');

// ---- memory: growth detaches typed arrays over the old buffer, and collected worldlines are released ----
const before = SRKit.wl.live(), pages = SRKit.raw.memory.buffer.byteLength;
let keep = [];
for (let i = 0; i < 3000; i++) { const s = []; for (let j = 0; j < 12; j++) s.push(mk.lin(j, j + 1, [i, j, 0], [0.1, 0, 0])); keep.push(SRKit.wl.make(s, 0, 0)); }
ok(SRKit.raw.memory.buffer.byteLength > pages, 'the module grew its memory');
near(S.at(kinds.hyp, 1).X[1], Math.SQRT2 - 1, 'and calls still read the right buffer after it grew');
ok(SRKit.wl.live() - before >= 3000, 'the module holds a copy of each worldline made');
keep = null;
if (typeof gc === 'function') {
  const settle = async () => { for (let i = 0; i < 20 && SRKit.wl.live() - before > 100; i++) { gc(); await new Promise(r => setTimeout(r, 20)); } };
  settle().then(() => { ok(SRKit.wl.live() - before <= 100, `collected worldlines are released (live ${SRKit.wl.live() - before})`); finish(); });
} else { console.log('     (run with --expose-gc to check that collected worldlines are released)'); finish(); }

function finish() {
  // ---- orbits and entries ----
  const entry = SRKit.plan.orbitEntry({ t: 0, tau: 0, X: [10, 0, 0] }, [0, 0, 0], 2, 0.6, g);
  ok(entry && entry.fromInside === false && entry.orb.r === 2 && Math.abs(entry.orb.w - 0.3) < 1e-15, 'orbit entry along a tangent');
  const circle = SRKit.plan.circPiece({ t: entry.end.t, tau: entry.end.tau, X: entry.end.X }, entry.orb, Infinity);
  const joined = SRKit.Worldline.from([...entry.pieces.map(p => p.g), circle.g], 0, 0);
  ok(joined.continuous(1e-9, 1e-9), 'the velocity is continuous where the run meets the circle'); joined.free();
  ok(SRKit.plan.orbitEntry({ t: 0, tau: 0, X: [0.1, 0, 0] }, [0, 0, 0], 2, 0.6, g).fromInside === true, 'orbit entry from inside');
  const exit = SRKit.plan.orbitExit({ t: 0, tau: 0, X: [2, 0, 0] }, [0, 0.6, 0], g);
  near(exit.dphi, Math.atanh(0.6), 'leaving an orbit spends its rapidity'); ok(exit.pieces[0].phase === 'brake' && SRKit.plan.brake(start, [1, 0, 0], 0, g).pieces.length === 0, 'brake');
  const plane = SRKit.plan.planeThrough([0, 0, 0], [3, 0, 4]); near(Math.hypot(...plane.e1), 1, 'planeThrough gives unit vectors'); near(plane.e1[0] * plane.e2[0] + plane.e1[1] * plane.e2[1] + plane.e1[2] * plane.e2[2], 0, 'at right angles', 1e-15);
  const o = SRKit.orbit.ofAcceleration(100, g); near(o.gamma, Math.hypot(1, Math.sqrt(g * 100)), 'loop gamma'); near(o.lapTau, o.period / o.gamma, 'lap proper time');
  near(SRKit.orbit.ofSpeed(100, o.beta).acc, g, 'the two ways to give an orbit agree', 1e-12);
  const os = SRKit.orbit.at(100, g, 0.3, 5); near(Math.hypot(...os.X), 100, 'on the circle'); near(SRKit.orbit.cone(100, g, 0.3, os.t, os.X, false), 5, 'an event on the orbit is received as it happens', 1e-9);

  // ---- rides ----
  const p = SRKit.trip.plan(1e5, g, 1e4, 1, true);
  ok(p.n === 4 && p.nOut === 2 && p.rt === true && p.leg.tc > 0 && Math.abs(p.gPeak - 1e4) < 1e-6, 'a ride with a cap, a stop and a return');
  const mid = SRKit.trip.at(p, p.tau * 0.3); near(SRKit.trip.tauU(p, mid.u), mid.tau, 'tauU inverts u', 1e-9); near(SRKit.trip.tauV(p, mid.v), mid.tau, 'tauV inverts v', 1e-6);
  ok(SRKit.trip.at(p, -1).t === -1 && SRKit.trip.at(p, p.tau + 2).t === p.tEnd + 2 && SRKit.trip.atClamped(p, p.tau + 2).t === p.tEnd, 'before launch, after arrival, and held within the ride');
  near(SRKit.trip.cone(p, mid.t, mid.x, 0, false), mid.tau, 'an event on the ride is received as it happens', 1e-9);
  assert.throws(() => SRKit.trip.at({ ...p }, 1), /SRKit\.trip\.plan/); checks++;
  for (const bad of [() => SRKit.trip.plan(-1, g), () => SRKit.trip.plan(1, 0), () => SRKit.trip.plan(1, g, 0.5), () => SRKit.trip.plan(1, g, Infinity, 1.5), () => SRKit.orbit.ofAcceleration(0, g)]) { assert.throws(bad, RangeError); checks++; }
  const fb = SRKit.trip.flipBurn(1e5, g); near(fb.tau, 22.37, 'oracle SR-03', 1e-3); near(fb.lag, 1.937, 'oracle: 1 g across the galaxy arrives 1.937 yr behind light', 1e-3);
  const r0 = SRKit.trip.reach(0, g, 10, 1.1); near(Math.hypot(...r0), 2 / g * (Math.cosh(g * 5) - 1), 'reach from rest is the flip-and-burn sphere', 1e-9);
  ok(SRKit.trip.reach(3, g, 1, 0) === null, 'reach is null when the budget cannot stop the ship');

  // ---- optics, rockets ----
  const fwd = SRKit.optics.observe(2, 1, 0, 0); near(fwd.D, Math.exp(2), 'dead ahead: Doppler e^phi'); near(SRKit.optics.observe(2, -1, 0, 0).D, Math.exp(-2), 'dead astern: e^-phi');
  const side = SRKit.optics.observe(2, 0, 1, 0); near(side.thp, 2 * Math.atan(Math.exp(-2)), 'aberration of a star abeam'); ok(SRKit.optics.observe(2, 0, 1, 0, false).thp === side.th, 'aberration can be switched off');
  const b = SRKit.optics.boost([0, 0, 1.5], 3, [1, 2, 0.5]), bb = SRKit.optics.boost([0, 0, -1.5], b.t, b.X); near(bb.t, 3, 'a boost and its inverse', 1e-13); near(bb.X[2], 0.5, 'a boost and its inverse (z)', 1e-13);
  near(SRKit.optics.compose([0.5, 0, 0], [0.5, 0, 0])[0], 0.8, 'velocities compose: 0.5 c and 0.5 c make 0.8 c');
  near(SRKit.optics.playback([0, 0, 0], [1, 0, 0], [0, 0, 0]), 1, 'two bodies at rest see each other in real time');
  const basis = SRKit.optics.basis([0.5, 0, 0]); ok(basis.f[0] === 0.5 && basis.u[2] === 1 && Math.abs(basis.r[1] + 0.5) < 1e-15, 'basis');
  near(SRKit.rocket.massRatio(2), Math.exp(2), 'photon rocket'); near(SRKit.rocket.phiOf(SRKit.rocket.massRatio(1.3, 0.1), 0.1), 1.3, 'the rocket equation inverts');
  near(SRKit.rocket.energy(1, 1, 0.5), 2 * Math.expm1(1), 'propellant energy'); near(SRKit.rocket.sailPhi(0.5), Math.log(2), 'sail');
  near(SRKit.kepler(1, 0), 1, 'Kepler on a circle'); near(SRKit.kepler(1, 0.5) - 0.5 * Math.sin(SRKit.kepler(1, 0.5)), 1, 'Kepler solves E - e sin E = M', 1e-14);
  near(SRKit.optics.resolution(1, 1e-6, 1e6), 1.22, 'Rayleigh'); near(SRKit.optics.luminosity(SRKit.units.R_SUN, 5772), SRKit.units.L_SUN, 'the Sun at 5,772 K', 2e-3);
  ok(SRKit.optics.bitRate(1e-9, 1e-6) > 5e9 && SRKit.optics.friis(1, 1, 1, 1, 1) > 0.6, 'link budget');

  // ---- cosmology and gravity, through the same module ----
  if (SRKit.flrw && SRKit.gr) {
    ok(Object.keys(SRKit.flrw).length > 5 && Object.keys(SRKit.gr).length > 5, 'flrw and gr are exposed');
    near(SRKit.flrw.age(), 13.79e9, 'the age of a Planck 2018 universe without radiation, in years', 5e-3);
    near(SRKit.gr.kerrIscoRadius(0), 6, 'Schwarzschild ISCO at 6 GM/c^2'); near(SRKit.gr.kerrIscoRadius(1), 1, 'extreme Kerr, prograde', 1e-6);
  }

  // ---- exact time: decimal strings or BigInt, the default clock ----
  const T = SRKit.time, epoch = '3' + '0'.repeat(110);
  ok(T.u384.bits === 384 && T.u512.bits === 512 && T.u768.bits === 768 && T.u1024.bits === 1024 && T.bits === 384 && T.ticksPerSecond === 1000 && T.gridPerMetre === 1, 'four widths; the default is 384 bits, milliseconds, metres');
  ok(T.u512.add(epoch, '1') === '3' + '0'.repeat(109) + '1' && T.u512.add(epoch, '-1') === '2' + '9'.repeat(110), 'add at an epoch no double can hold');
  ok(T.add(10n ** 100n, 5n) === 10n ** 100n + 5n && T.diff(7n, 10n ** 30n) === 7n - 10n ** 30n, 'BigInt in, BigInt out');
  ok(T.cmp('100', '99') === 1 && T.cmp('5', '5') === 0 && T.cmp(3n, 10n ** 90n) === -1, 'cmp');
  ok(T.mod('1' + '0'.repeat(100), 31557600000) === Number(10n ** 100n % 31557600000n), 'mod is the exact phase');
  ok(T.seconds('1000000000000000001500', '1000000000000000000000') === 1.5 && T.seconds('1' + '0'.repeat(40), '0') === null, 'seconds: exact subtraction first, null when too far apart');
  ok(T.lightTicks([299792458, 0, 0]) === '1000' && T.lightTicks([299792459, 0, 0]) === '1001' && T.lightTicks([3 * 299792458, 4 * 299792458, 0]) === '5000', 'light ticks round up, exactly');
  const arrival = T.u512.firstArrival(epoch, [4.0175e16, 0, 0]);
  ok(T.u512.arrived(epoch, arrival, [4.0175e16, 0, 0]) && !T.u512.arrived(epoch, T.u512.add(arrival, '-1'), [4.0175e16, 0, 0]), 'arrived flips exactly at the first arrival');
  ok(T.classify('1000', [299792458, 0, 0]) === 'null' && T.classify('1000', [299792457, 0, 0]) === 'timelike' && T.classify('-1000', [0, 299792459, 0]) === 'spacelike', 'classify, exactly on the cone');
  assert.throws(() => T.add('-5', '1'), RangeError); assert.throws(() => T.add('1' + '0'.repeat(120), '1'), /out of range/); assert.throws(() => T.u512.add('0', '-1'), /out of range/); checks += 3;

  // ---- any clock: four widths, each at any tick and any grid (R71.4) ----
  const C = 299792458n, RES = [[1e3, 1], [1e6, 1], [1e9, 1e3], [1e12, 1e6], [1, 1], [60, 100], [1024, 3]];
  for (const bits of [384, 512, 768, 1024]) {
    for (const [tps, gpm] of RES) {
      const K = T.clock({ bits, ticksPerSecond: tps, gridPerMetre: gpm }), one = [C * BigInt(gpm), 0n, 0n], sent = 3n * 10n ** 60n;
      const heard = K.firstArrival(sent, one);
      ok(heard - sent === BigInt(tps) && K.lightTicks(one) === String(tps), `U${bits} at ${tps} ticks/s on a 1/${gpm} m grid: one light-second is one second of ticks`);
      ok(K.arrived(sent, heard, one) && !K.arrived(sent, heard - 1n, one) && K.classify(BigInt(tps), one) === 'null' && K.classify(BigInt(tps) - 1n, one) === 'spacelike', `U${bits} ${tps}/${gpm}: arrived flips at that tick, which is on the cone`);
      ok(K.lightTicks([C * BigInt(gpm) + 1n, 0n, 0n]) === String(tps + 1), `U${bits} ${tps}/${gpm}: one grid step further is one tick later`);
      ok(K.lightTicks([299792458, 0, 0]) === String(tps), `U${bits} ${tps}/${gpm}: the same from metres`);
    }
  }
  const ns768 = T.clock({ bits: 768, ticksPerSecond: 1e9, gridPerMetre: 1e3 });
  near(ns768.reachLog10Years, 214.7, 'U768 at a nanosecond reaches 10^214.7 years', 1e-3); near(T.reachLog10Years, 105.1, 'U384 at a millisecond reaches 10^105.1 years', 1e-3);
  near(ns768.lightStepsPerTick, 299.792458, 'light covers 300 mm a nanosecond', 1e-12); near(T.clock({ ticksPerSecond: 1e9 }).lightStepsPerTick, 0.299792458, 'and less than one metre: a nanosecond needs a finer grid', 1e-12);
  ok(T.u1024.max === String(2n ** 1024n - 1n) && T.u768.max === String(2n ** 768n - 1n) && T.max === String(2n ** 384n - 1n), 'the largest time of each width, 309 digits at 1,024 bits');
  ok(T.u1024.add(T.u1024.max, '-1') === String(2n ** 1024n - 2n) && T.u1024.diff('0', T.u1024.max) === '-' + T.u1024.max, 'arithmetic at the full width');
  assert.throws(() => T.u768.add(T.u768.max, '1'), /out of range/); assert.throws(() => T.u768.add('1' + '0'.repeat(232), '0'), /out of range/); checks += 2;
  // 1e200 years of picoseconds: beyond 512 bits, exact in 768
  const ps = T.clock({ bits: 768, ticksPerSecond: 1e12, gridPerMetre: 1e6 }), far = 31557600n * 10n ** 212n;
  assert.throws(() => T.clock({ bits: 512, ticksPerSecond: 1e12, gridPerMetre: 1e6 }).add(far, 1n), /out of range/); checks++;
  const across = [88n * 10n ** 29n, 3n, -7n], echo = ps.firstArrival(far, across);
  ok(ps.arrived(far, echo, across) && !ps.arrived(far, echo - 1n, across) && typeof echo === 'bigint', 'an arrival across 8.8e24 m in micrometres, exact to the picosecond at 1e200 years');
  // far positions as grid steps are exact; as floats they are not: 4.0175e16 m + 1 m is not a double
  ok(T.lightTicks([40175000000000001n, 0n, 0n]) !== T.lightTicks([40175000000000001n + 299792n, 0n, 0n]) && T.lightTicks(['40175000000000000', '0', '0']) === T.lightTicks([4.0175e16, 0, 0]), 'a displacement in grid steps, as BigInt or text');
  assert.throws(() => T.lightTicks([1, 2n, 3]), TypeError); assert.throws(() => T.lightTicks(['1.5', '0', '0']), RangeError); assert.throws(() => T.lightTicks([10n ** 40n, 0n, 0n]), /out of range/); checks += 3;
  // seconds and phase follow the tick
  ok(ns768.seconds(10n ** 150n + 90n * 10n ** 9n, 10n ** 150n) === 90 && ns768.seconds('0', '1500000000') === -1.5 && ns768.seconds(2n ** 60n, 0n) === null, 'seconds at a nanosecond tick');
  const yearNs = 31557600n * 10n ** 9n;   // 3.2e16 ticks: past 2^53, so the period is a BigInt
  ok(ns768.mod(10n ** 150n, yearNs) === 10n ** 150n % yearNs && ns768.mod(10n ** 150n, 1000) === 0 && ns768.mod('7', 5) === 2 && T.mod('7', 0) === null && T.mod('7', 2.5) === null && T.mod('7', 2 ** 60) === null, 'mod with a period past 2^53 as a BigInt');
  const yearPs = 31557600n * 10n ** 12n;   // a year of picoseconds does not fit 64 bits
  ok(ps.mod(far + 12345n, yearPs) === 12345n && ps.mod(far, 0n) === null && ps.mod(far, -5n) === null, 'mod with a period past 2^64');
  for (const bad of [{ bits: 256 }, { bits: 640 }, { ticksPerSecond: 0 }, { ticksPerSecond: 1e13 }, { ticksPerSecond: 1000.5 }, { gridPerMetre: 0 }, { gridPerMetre: 1e7 }, { ticksPerSecond: NaN }]) { assert.throws(() => T.clock(bad), RangeError); checks++; }
  ok(Object.isFrozen(ns768) && T.clock().bits === 384 && T.clock().ticksPerSecond === 1000, 'clock() with no options is the default clock');

  // ---- other hosts ----
  const src = fs.readFileSync(path.join(dist, 'srkit.js'), 'utf8');
  const bare = vm.createContext({}); vm.runInContext(src, bare);   // no console, no Buffer, no atob, no require
  ok(vm.runInContext('SRKit.seg.at(SRKit.seg.mk.hyp(0, 1, [0, 0, 0], [1, 0, 0], 1, 0), 1).X[0]', bare) === S.at(mk.hyp(0, 1, [0, 0, 0], [1, 0, 0], 1, 0), 1).X[0], 'runs in a bare vm context, with the same bits');
  assert.throws(() => vm.runInContext(src, vm.createContext({ WebAssembly: undefined })), /WebAssembly/); checks++;
  const second = SRKit.create(wasm);
  ok(second !== SRKit && second.info.bytes === wasm.length && second.trip.flipBurn(1, 1).tau === SRKit.trip.flipBurn(1, 1).tau && second.calls() < SRKit.calls(), 'create(bytes) makes an independent instance');
  ok(SRKit.calls() > 5000, 'the module did the work');
  console.log(`srkit.js: ${checks} checks passed (module ${SRKit.info.bytes} bytes, ${SRKit.calls()} calls, node ${process.version})`);
}
