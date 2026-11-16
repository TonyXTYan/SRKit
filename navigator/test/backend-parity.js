'use strict';
// Run from any directory: node navigator/test/backend-parity.js
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
function fields(actual, expected, names, label, rel, abs) {
  for (const name of names) near(actual[name], expected[name], `${label}.${name}`, rel, abs);
}
function vector(actual, expected, label, rel, abs) {
  assert.equal(actual.length, expected.length, label);
  actual.forEach((value,i)=>near(value,expected[i],`${label}[${i}]`,rel,abs));
}
for (const call of [()=>N.trip(-1,1,Infinity,0,false),()=>N.trip(1,0,Infinity,0,false),()=>N.trip(1,1,.5,0,false),()=>N.orbit(0,1)]) {
  assert.throws(call, /SRKit|invalid|range|input/i, 'invalid physical inputs must fail explicitly');
  checks++;
}
const startCalls = N.calls();
let plans = 0;
for (const d of [4.2465, 1e5, 2.5e6]) for (const g of [.01, 1]) {
  for (const cap of [Infinity, 1e5]) for (const stops of [0, 2]) for (const rt of [false, true]) {
    const p = N.trip(d,g,cap,stops,rt), q = legacy.trip(d,g,cap,stops,rt);
    const label = `trip(${[d,g,cap,stops,rt].join(',')})`;
    fields(p,q,['d','a','n','nOut','legTau','tau','du','dv','uEnd','vEnd','tEnd','tOut','xEnd','lag','gPeak','phc'],label,3e-10,1e-8);
    fields(p.leg,q.leg,['a','phc','ta','tc'],`${label}.leg`,3e-10,1e-8);
    assert.equal(p.rt,q.rt);
    for (const fraction of [-.2,0,.02,.2,.5,.75,.98,1,1.2]) {
      const tau = q.tau*fraction;
      const actual=N.at(p,tau),expected=legacy.at(q,tau);
      fields(actual,expected,['tau','u','v','t','x','phi','gamma','Dh','Dd'],`${label}@${fraction}`,2e-9,2e-8);
      // A one-ulp plan-duration difference can place an exact waypoint on
      // either side of a rest-to-rest join. The heading is arbitrary at rest;
      // require direction parity whenever either ship is actually moving.
      if(Math.max(Math.abs(actual.phi),Math.abs(expected.phi))>1e-10) near(actual.dir,expected.dir,`${label}@${fraction}.dir`);
      near(actual.dir*Math.tanh(actual.phi),expected.dir*Math.tanh(expected.phi),`${label}@${fraction}.velocity`,2e-9,2e-8);
    }
    for (const fraction of [-.1,0,.1,.35,.7,.99,1,1.1]) {
      near(N.tauU(p,q.uEnd*fraction),legacy.s1TauU(q,q.uEnd*fraction),`${label}.tauU(${fraction})`,2e-8,2e-7);
      near(N.tauV(p,q.vEnd*fraction),legacy.s1TauV(q,q.vEnd*fraction),`${label}.tauV(${fraction})`,2e-8,2e-7);
    }
    // Off-axis photons test the actual receiver/ping root solver as well as
    // the analytic inverse used by the on-axis charts.
    if (stops === 0 && cap === Infinity) {
      for (const t of [-2, q.tEnd*.3, q.tEnd+5]) for (const send of [false,true]) {
        const along=.3*d, perp=.1*d;
        near(N.recv(p,t,along,perp,send),legacy.recv(q,t,along,perp,send),`${label}.recv(${t},${send})`,2e-8,2e-7);
      }
    }
    plans++;
  }
}
for (const r of [1,100,50000]) for (const g of [.01,1]) {
  const L=N.orbit(r,g), R=legacy.orbit(r,g);
  fields(L,R,['a','gamma','beta','phi','lapT','lapTau','omega','sagnac'],`orbit(${r},${g})`);
  for (const fraction of [-.25,0,.1,.5,1,1.25]) {
    const a0=.37,tau=R.lapTau*fraction;
    const a=N.orbitAt(L,a0,tau),b=legacy.orbitAt(R,a0,tau);
    fields(a,b,['t','al'],'orbit state',3e-11,1e-9);
    vector(a.X,b.X,'orbit position',3e-11,1e-9);
    vector(a.V,b.V,'orbit velocity');
    vector(a.heading,b.heading,'orbit heading');
  }
  for (const P of [[0,0,0],[1.2*r,-.5*r,.2*r]]) for (const t of [0,R.lapT*.2]) {
    for (const send of [false,true]) near(N.orbitRecv(L,.37,t,P,send),legacy.orbitRecv(R,.37,t,P,send),'orbit light-cone intersection',5e-10,1e-8);
  }
}
for (const phi of [0,.5,5,14]) for (const theta of [0,.1,Math.PI/2,Math.PI-.1,Math.PI]) {
  for (const aberr of [false,true]) {
    const a=Math.cos(theta),br=Math.sin(theta)*.6,bu=Math.sin(theta)*.8;
    const expected=legacy.optics(phi,a,br,bu,aberr),actual=N.observe(phi,a,br,bu,aberr);
    fields(actual,expected,['D','th','thp','f','r','u'],`optics(${phi},${theta},${aberr})`,3e-10,1e-12);
    near(actual.f**2+actual.r**2+actual.u**2,1,'observed unit direction',1e-10,1e-12);
  }
}
for (const w of [[0,0,0],[.5,0,0],[.4,-.3,.2],[5,-1,3]]) {
  const dt=3.7,X=[2,-.2,.9], actual=N.boost(w,dt,X),expected=legacy.boost(w,dt,X);
  near(actual.t,expected.t,'boost time',1e-10,1e-10);
  vector(actual.X,expected.X,'boost position',1e-10,1e-10);
  const inverse=N.boost(w.map(x=>-x),actual.t,actual.X);
  near(inverse.t,dt,'inverse boost time',1e-8,1e-8);
  vector(inverse.X,X,'inverse boost position',1e-8,1e-8);
}

assert.ok(N.calls()>startCalls+100,'parity suite must invoke compiled WebAssembly');
// An absent runtime must fail clearly, never switch to the legacy JS fixture.
const missing = vm.createContext({console, WebAssembly: undefined});
assert.throws(()=>vm.runInContext(source, missing), /WebAssembly|wasm|SRKit/i);
checks++;
console.log(`SRKit Navigator parity: ${checks} assertions passed across ${plans} trip plans; ${N.calls()-startCalls} Wasm calls`);
