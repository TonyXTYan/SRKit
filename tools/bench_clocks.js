// U384 against U512 as the stored clock, in WebAssembly (question R23.8). Run after tools/build_js.py:
//     node srkit/tools/bench_clocks.js
// Each round advances a time by one tick, takes a difference, compares, reduces to an orbital phase and makes an
// exact arrival check. The loop runs inside the module, so this measures the arithmetic, not the JS boundary.
'use strict';
const SRKit = require('../dist/srkit.js');
const n = 300000, rounds = 7;
for (const width of [384, 512]) {
  const times = [];
  let checksum;
  for (let r = 0; r < rounds; r++) {
    const t0 = process.hrtime.bigint();
    checksum = SRKit.raw.srkit_bench_clock(width, n);
    times.push(Number(process.hrtime.bigint() - t0) / n);
  }
  times.sort((a, b) => a - b);
  console.log(`U${width}: ${times[rounds >> 1].toFixed(1)} ns per round (median of ${rounds}; best ${times[0].toFixed(1)}), checksum ${checksum}`);
}
console.log(`module ${SRKit.info.bytes} bytes, sha256 ${SRKit.info.sha256.slice(0, 16)}…; node ${process.version}`);
