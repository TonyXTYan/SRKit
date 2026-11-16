// The stored clock at three widths and three tick lengths, in WebAssembly (questions R23.8 and R71.4). Run after
// tools/build_js.py:
//     node srkit/tools/bench_clocks.js
// Each round advances a time by one tick, takes a difference, compares, reduces to a yearly phase and makes an
// exact arrival check. The loop runs inside the module, so this measures the arithmetic, not the JS boundary.
'use strict';
const SRKit = require('../dist/srkit.js');
const n = 200000, rounds = 7;
console.log('ns per round in WebAssembly (median of ' + rounds + ' batches of ' + n + '; best in brackets)');
for (const [tick, tps, gpm, grid] of [['ms', 1e3, 1, 'metre'], ['us', 1e6, 1, 'metre'], ['ns', 1e9, 1e3, 'millimetre']]) {
  for (const width of [384, 512, 1024]) {
    const times = [];
    let checksum;
    for (let r = 0; r < rounds; r++) {
      const t0 = process.hrtime.bigint();
      checksum = SRKit.raw.srkit_bench_clock(width, tps, gpm, n);
      times.push(Number(process.hrtime.bigint() - t0) / n);
    }
    times.sort((a, b) => a - b);
    console.log(`U${String(width).padEnd(5)}${tick}  ${grid.padEnd(11)} ${times[rounds >> 1].toFixed(0).padStart(6)}  (${times[0].toFixed(0)})   checksum ${checksum}`);
  }
}
console.log(`module ${SRKit.info.bytes} bytes, sha256 ${SRKit.info.sha256.slice(0, 16)}…; node ${process.version}`);
