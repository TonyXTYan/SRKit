'use strict';
// Run from any directory: node navigator/test/reach-traveller.js
// This compares the compiled WebAssembly path with a preserved independent JS
// reference. It intentionally never supplies a production physics fallback.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const legacy = require('./fixtures/navigator-legacy.js');
const backendPath = path.resolve(__dirname, '../src/srkit-backend.js');
const source = fs.readFileSync(backendPath, 'utf8');
const context = vm.createContext({console});
vm.runInContext(source, context, {filename: backendPath});
const N = context.SRKitNavigator;
assert.ok(N && N.info, 'SRKit backend must initialize synchronously');
assert.ok(/srkit/i.test(JSON.stringify(N.info)), 'backend identifies SRKit');
assert.equal(typeof N.calls, 'function', 'backend exposes a Wasm call counter');
let checks = 0;
function near(actual, expected, label, rel = 3e-11, abs = 2e-12) {
  assert.ok(Number.isFinite(actual) && Number.isFinite(expected), `${label}: finite values required (${actual}, ${expected})`);
  assert.ok(Math.abs(actual-expected) <= abs+rel*Math.max(Math.abs(actual),Math.abs(expected)), `${label}: ${actual} versus ${expected}`);
  checks++;
}
function vector(actual, expected, label, rel, abs) {
  assert.equal(actual.length, expected.length, label);
  actual.forEach((value,i)=>near(value,expected[i],`${label}[${i}]`,rel,abs));
}
const startCalls = N.calls();
// The reach contour must agree with the original go/flip/brake construction,
// and its rest-frame limit is a symmetric flip-and-burn sphere.
for (const phi of [0,.5,3,8]) for (const a of [.01,legacy.G_LY]) {
  const stop=phi/a;
  if(phi>0) {
    assert.equal(N.reach(phi,a,stop*.75,.3),null,'insufficient budget cannot stop the ship');
    checks++;
  }
  for(const budget of [stop,stop+.1,stop+2,stop+20]) for(const heading of [-Math.PI,-Math.PI/2,0,.2,Math.PI/2,Math.PI]) {
    const R={phi,a,g0:Math.cosh(phi),gb0:Math.sinh(phi)};
    const expected=legacy.reachAS(R,budget,heading),actual=N.reach(phi,a,budget,heading);
    assert.ok(expected && actual,'valid reach budget returns a point');
    vector(actual,expected,`reach(${phi},${a},${budget},${heading})`,3e-9,1e-7);
  }
}
for(const a of [.01,legacy.G_LY]) for(const budget of [0,.1,2,20]) for(const heading of [-2,-.4,0,1,Math.PI]) {
  const point=N.reach(0,a,budget,heading);
  const radius=4/a*Math.sinh(a*budget/4)**2;
  vector(point,[radius*Math.cos(heading),radius*Math.sin(heading)],'rest reach equals flip-and-burn radius',3e-10,1e-9);
}
// Independent monotone light-cone bisection for a traveller waiting until t=0
// then coasting in a full 3D direction. It does not share the quadratic solve.
function travellerReference(beta,heading,observerT,observer) {
  const length=Math.hypot(...heading),direction=heading.map(v=>v/length);
  const residual=te=>observerT-te-Math.hypot(...observer.map((v,i)=>v-direction[i]*beta*Math.max(0,te)));
  let hi=observerT,step=Math.max(1,Math.abs(observerT),Math.hypot(...observer)),lo=hi-step;
  for(let i=0;residual(lo)<0 && i<100;i++) {step*=2;lo=hi-step;}
  assert.ok(residual(lo)>=0,'independent traveller reference has a bracket');
  for(let i=0;i<200;i++) {const mid=lo+(hi-lo)/2;if(mid===lo||mid===hi)break;if(residual(mid)>=0)lo=mid;else hi=mid;}
  return lo+(hi-lo)/2;
}
for(const beta of [0,.5,.9,.99]) for(const heading of [[1,0,0],[1,2,-3],[-.3,.1,.7]]) {
  for(const observerT of [-20,0,2,100,1e5]) for(const observer of [[0,0,0],[3,-4,5],[-7,2,1]]) {
    const actual=N.travellerEmission(beta,heading,observerT,observer);
    const expected=travellerReference(beta,heading,observerT,observer);
    near(actual,expected,'traveller emission matches independent 3D light-cone bisection',2e-10,2e-9);
    assert.ok(actual<=observerT+1e-9,'traveller emission cannot follow observation');
    checks++;
    if(observerT-Math.hypot(...observer)<=0) near(actual,observerT-Math.hypot(...observer),'traveller waits at rest before launch',1e-12,1e-10);
  }
}

assert.ok(N.calls()>startCalls+100,'reach/traveller checks must execute compiled WebAssembly');
console.log(`SRKit Navigator reach/traveller parity: ${checks} assertions passed; ${N.calls()-startCalls} Wasm calls`);
