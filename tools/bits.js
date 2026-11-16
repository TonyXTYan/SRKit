// The same results, to the bit, from the native build and the WebAssembly build.
//     python3 -I tools/build_js.py && node tools/bits.js
// Runs `cargo run -p srkit-wasm --example bits` (native: whatever this machine is) and the identical list of calls
// through dist/srkit.wasm, and compares bit patterns. Every transcendental goes through the pure-Rust libm, so the two
// must agree exactly. Run it on each architecture that matters; a difference is a determinism bug.
'use strict';
const { execFileSync } = require('node:child_process'), path = require('node:path');
const x = require('../dist/srkit.js').raw, root = path.resolve(__dirname, '..');
const native = execFileSync('cargo', ['run', '-q', '--release', '--locked', '-p', 'srkit-wasm', '--example', 'bits'], { cwd: root, encoding: 'utf8' }).trim().split('\n');
const seg = (v) => { for (let i = 0; i < 21; i++) x.srkit_set(i, v[i] === undefined ? 0 : v[i]); };
const get = (a, b) => { const o = []; for (let i = a; i < b; i++) o.push(x.srkit_get(i)); return o; };
const g = x.unit_constant(14), out = [g], INF = Infinity;
for (const k of [[2, 2000, INF, 0.1, -0.2, 0.3, 0.6, 0, 0.8, g, 0], [2, 2000, 2010, 0.1, -0.2, 0.3, 0.6, 0, 0.8, -0.37, 4.5], [3, 2000, INF, 1, 2, 3, 7.5, 0.11, 1, 0, 0, 0, 1, 0, 0.7],
  [5, 2000, INF, 0, 0, 0, 1.2, -0.4, 0.3, 0.2, 0.9, -0.1, 0.83], [6, 2000, INF, 0, 0, 0, 0, 1, 0, 0.3, 0.71, 0.4]]) {
  seg(k);
  for (const t of [2000.001, 2001.7, 2009.3]) { out.push(x.seg_dtau(t)); x.seg_at(t); out.push(...get(0, 6)); out.push(x.seg_t_of_tau(x.seg_dtau(t))); }
  out.push(x.seg_recv(1999, 3, -2, 1, 2000, INF), x.seg_send(2030, 3, -2, 1, -INF, INF));
}
for (const [d, cap, stops, rt] of [[4.2465, INF, 0, 0], [1e5, 1e4, 1, 1], [2.5e6, INF, 0, 0]]) {
  x.trip_plan(d, g, cap, stops, rt); out.push(...get(5, 21));
  const tau = x.srkit_get(6) * 0.37;
  x.trip_state(d, g, cap, stops, rt, tau, 0); out.push(...get(1, 10));
  out.push(x.trip_tau_of(d, g, cap, stops, rt, d * 0.5, 1), x.trip_light_cone(d, g, cap, stops, rt, d * 0.4, d * 0.3, d * 0.01, 0));
}
x.reach_point(1.3, g, 9, 0.8); out.push(...get(0, 2));
x.observe(3.3, 0.2, -0.5, 0.7, 1); out.push(...get(0, 6));
x.boost(0.3, -1.1, 2, 5, 1, 2, 3); out.push(...get(0, 4));
out.push(x.rocket_mass_ratio(3.7, 0.3), x.sail_rapidity(0.42), x.kepler_eccentric_anomaly(2.2, 0.6), x.black_body_luminosity(6.957e8, 5772));
const h0 = x.flrw_planck_2018(0), om = x.flrw_planck_2018(1);
for (const [code, a, b] of [[0, 0, 0], [1, 5e9, 0], [4, 13.8e9, 0], [5, 0, 0], [6, 13.8e9, 0], [11, 1e10, 0], [13, 13.8e9, 5e9], [20, g, 0], [23, 0, 0]]) out.push(x.flrw(code, h0, om, a, b, 0));
for (const [code, a] of [[0, 0], [0, 0.9], [0, -0.5], [1, 0.9], [2, 0.9], [3, 1e9], [4, 1e9], [12, 0.7]]) out.push(x.gr(code, a, 0, 0, 0));
out.push(x.gr(9, 1.98841e30, 1.496e11, 2.279e11, 6.957e8));
const hex = (v) => { const b = new DataView(new ArrayBuffer(8)); b.setFloat64(0, v); return b.getBigUint64(0).toString(16).padStart(16, '0'); };
let bad = 0;
if (native.length !== out.length) { console.log(`FAIL the two lists differ in length: native ${native.length}, wasm ${out.length}`); bad++; }
for (let i = 0; i < Math.min(native.length, out.length); i++) if (native[i] !== hex(out[i])) { if (bad < 10) console.log(`FAIL result ${i}: native ${native[i]}, wasm ${hex(out[i])} (${out[i]})`); bad++; }
console.log(bad ? `${bad} of ${out.length} results differ` : `${out.length} results identical to the bit: native ${process.arch} and wasm32 (node ${process.version})`);
process.exitCode = bad ? 1 : 0;
