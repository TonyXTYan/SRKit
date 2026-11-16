'use strict';
// Independent reference: unchanged physics from claude/maps-and-navigation
// commit e03e19c of the game's repository (docs/design/mockups/src/navigator.src.html there).
// Retained only for regression tests; production uses SRKit WebAssembly.
const G_LY = 9.80665 * (365.25 * 86400) / 299792458;
const dot = (a,b) => a.reduce((v,x,i) => v+x*b[i],0);
const len = a => Math.hypot(...a);
const sub = (a,b) => a.map((x,i)=>x-b[i]);
const scl = (a,s) => a.map(x=>x*s);

  function legPlan(L, a, cap) {
    const gPk = 1 + a * L / 2;
    if (cap < gPk) { const phc = Math.acosh(cap), dacc = (cap - 1) / a; return { a, phc, ta: phc / a, tc: (L - 2 * dacc) / Math.sinh(phc) }; }
    const ph = Math.acosh(gPk); return { a, phc: ph, ta: ph / a, tc: 0 };
  }
  // a ride of nOut identical flip-and-burn legs out to distance d; a round trip (rt) adds the same legs back to the start.
  // u = t - x and v = t + x are the null coordinates: the home time and the destination time (+ d) being seen.
  function trip(d, g, cap, stops, rt) {
    const a = g * G_LY, nOut = stops + 1, n = rt ? 2 * nOut : nOut, leg = legPlan(d / nOut, a, cap);
    const ea = Math.exp(leg.phc), eb = Math.exp(-leg.phc);
    const du = (2 / a) * (1 - eb) + eb * leg.tc, dv = (2 / a) * (ea - 1) + ea * leg.tc, legTau = 2 * leg.ta + leg.tc;
    const tEnd = n * (du + dv) / 2;
    return { d, a, n, nOut, rt: !!rt, leg, legTau, tau: n * legTau, du, dv, uEnd: rt ? tEnd : n * du, vEnd: rt ? tEnd : n * dv, tEnd,
      tOut: nOut * (du + dv) / 2, xEnd: rt ? 0 : d, lag: nOut * (du + dv) / 2 - d, gPeak: Math.cosh(leg.phc), phc: leg.phc };
  }
  // Ship time may run before launch (you wait at home, at rest: ship time is galactic time) and after the ride (at rest where it
  // ends: both clocks tick together from the arrival event)
  function at(p, tau) {
    if (tau < 0) return { tau, u: tau, v: tau, t: tau, x: 0, phi: 0, dir: 1, gamma: 1, Dh: 1, Dd: 1 };
    if (tau > p.tau) { const e = atRide(p, p.tau), dt = tau - p.tau; return { ...e, tau, u: e.u + dt, v: e.v + dt, t: e.t + dt }; }
    return atRide(p, tau);
  }
  function atRide(p, tau) {
    tau = Math.max(0, Math.min(p.tau, tau));
    if (p.tau - tau <= 1e-9 * p.tau) tau = p.tau;   // a slider or float a hair short of the end is arrival (exactly at rest)
    const i = Math.min(p.n - 1, Math.floor(tau / p.legTau)), s = tau - i * p.legTau, { a, phc, ta, tc } = p.leg;
    const ea = Math.exp(phc), eb = Math.exp(-phc), back = p.rt && i >= p.nOut, j = back ? i - p.nOut : i;
    let lu = 0, lv = 0, phi;   // this leg's own null-coordinate increments (along its direction of travel)
    if (s <= ta) { phi = a * s; lu = (1 - Math.exp(-phi)) / a; lv = (Math.exp(phi) - 1) / a; }
    else if (s <= ta + tc) { phi = phc; lu = (1 - eb) / a + eb * (s - ta); lv = (ea - 1) / a + ea * (s - ta); }
    else { const s2 = Math.min(s - ta - tc, ta); phi = tau >= p.tau || s2 >= ta ? 0 : phc - a * s2;
      lu = (1 - eb) / a + eb * tc + (Math.exp(-phi) - eb) / a; lv = (ea - 1) / a + ea * tc + (ea - Math.exp(phi)) / a; }
    // heading home, the roles swap: home's light now comes from ahead (u grows like v did), the destination's from behind
    const u = back ? p.nOut * p.du + j * p.dv + lv : i * p.du + lu, v = back ? p.nOut * p.dv + j * p.du + lu : i * p.dv + lv;
    const ed = Math.exp(phi), em = Math.exp(-phi);
    return { tau, u, v, t: (u + v) / 2, x: (v - u) / 2, phi, dir: back ? -1 : 1, gamma: Math.cosh(phi), Dh: back ? ed : em, Dd: back ? em : ed };
  }

  function s1Ship(p, tau) {   // { u, v, ph }: ph your rapidity along +x (negative heading home); at rest at home before launch and after arrival
    if (tau <= 0) return { u: tau, v: tau, ph: 0 };
    if (tau >= p.tau) { const e = tau - p.tau; return { u: p.uEnd + e, v: p.vEnd + e, ph: 0 }; }
    const s = at(p, tau); return { u: s.u, v: s.v, ph: s.dir * s.phi };
  }
  function s1Leg(p, l, k) {   // ship time into one leg at which its own null increment (k 'u' or 'v', along its travel) reaches l
    const { a, phc, ta, tc } = p.leg, ea = Math.exp(phc), eb = Math.exp(-phc);
    if (k === 'u') { const A = (1 - eb) / a, B = A + eb * tc;
      if (l <= A) return -Math.log1p(-Math.min(a * Math.max(0, l), 1 - eb)) / a;
      if (l <= B) return ta + (l - A) / eb;
      return ta + tc + (phc + Math.log(Math.min(1, eb + a * (l - B)))) / a; }
    const A = (ea - 1) / a, B = A + ea * tc;
    if (l <= A) return Math.log1p(a * Math.max(0, l)) / a;
    if (l <= B) return ta + (l - A) / ea;
    return ta + tc + (phc - Math.log(Math.max(1, ea - a * (l - B)))) / a;
  }
  function s1TauU(p, U) {   // the ship time at which your u = t - x reaches U (a light ray heading for the destination meets you)
    if (U <= 0) return U;
    if (U >= p.uEnd) return p.tau + (U - p.uEnd);
    const n = p.nOut;
    if (!p.rt || U < n * p.du) { const i = Math.min(n - 1, Math.floor(U / p.du)); return i * p.legTau + s1Leg(p, U - i * p.du, 'u'); }
    const W = U - n * p.du, j = Math.min(n - 1, Math.floor(W / p.dv)); return (n + j) * p.legTau + s1Leg(p, W - j * p.dv, 'v');
  }
  function s1TauV(p, V) {   // the ship time at which your v = t + x reaches V (a light ray heading home meets you)
    if (V <= 0) return V;
    if (V >= p.vEnd) return p.tau + (V - p.vEnd);
    const n = p.nOut;
    if (!p.rt || V < n * p.dv) { const i = Math.min(n - 1, Math.floor(V / p.dv)); return i * p.legTau + s1Leg(p, V - i * p.dv, 'v'); }
    const W = V - n * p.dv, j = Math.min(n - 1, Math.floor(W / p.du)); return (n + j) * p.legTau + s1Leg(p, W - j * p.du, 'u');
  }

  function s1Newton(f, lo, hi, x0, tol = 1e-15) {   // f(x) -> [value, slope], increasing, f(lo) <= 0 <= f(hi): Newton, kept inside the bracket
    let a = lo, b = hi, x = Math.min(Math.max(x0, a), b);
    for (let i = 0; i < 60; i++) {
      const [fv, fp] = f(x); if (fv === 0) return x; if (fv > 0) b = x; else a = x;
      let xn = x - fv / fp; if (!(xn > a && xn < b)) xn = 0.5 * (a + b);
      if (Math.abs(xn - x) <= tol * Math.max(1, Math.abs(x)) || b - a <= 2 * tol * Math.max(1, Math.abs(x))) return xn;
      x = xn;
    }
    return x;
  }
  // The ship time at which the light of event (t, P) reaches you. With a = P.m and b its distance from the line, it lies between
  // the receptions of the line's events (t, a) and (t + b, a), which are closed form; between them g = t_ship - t - r is increasing.
  function s1Recv(G, t, P) {
    const p = G.p, a = dot(P, G.m), b = len(sub(P, scl(G.m, a))), line = (tt) => Math.max(s1TauU(p, tt - a), s1TauV(p, tt + a)), lo = line(t);
    if (!(b > 1e-13 * (1 + Math.abs(a)))) return lo;
    return s1Newton((tau) => { const s = s1Ship(p, tau), xs = (s.v - s.u) / 2, r = Math.hypot(a - xs, b); return [(s.u + s.v) / 2 - t - r, Math.cosh(s.ph) - (xs - a) * Math.sinh(s.ph) / r]; }, lo, line(t + b), lo);
  }
  // The ship time of the ping that reaches event (t, P): between the pings of the line's events (t - b, a) and (t, a)
  function s1Ping(G, t, P) {
    const p = G.p, a = dot(P, G.m), b = len(sub(P, scl(G.m, a))), line = (tt) => Math.min(s1TauU(p, tt - a), s1TauV(p, tt + a)), hi = line(t);
    if (!(b > 1e-13 * (1 + Math.abs(a)))) return hi;
    return s1Newton((tau) => { const s = s1Ship(p, tau), xs = (s.v - s.u) / 2, r = Math.hypot(a - xs, b); return [(s.u + s.v) / 2 + r - t, Math.cosh(s.ph) + (xs - a) * Math.sinh(s.ph) / r]; }, line(t - b), hi, hi);
  }


function orbit(r,g) {
  const a=g*G_LY,gb=Math.sqrt(a*r),gamma=Math.sqrt(1+gb*gb),beta=gb/gamma;
  const lapT=2*Math.PI*r/beta;
  return {r,a,gamma,beta,phi:Math.acosh(gamma),lapT,lapTau:lapT/gamma,omega:beta/r,sagnac:2*Math.PI*r*beta*gamma};
}
function orbitAt(L,a0,tau) {
  const t=L.gamma*tau,al=a0+L.omega*t,c=Math.cos(al),s=Math.sin(al);
  return {t,X:[L.r*c,L.r*s,0],V:[-L.beta*s,L.beta*c,0],heading:[-s,c,0],al};
}
function optics(phi,a,br,bu,aberr=true) {
  const transverse=Math.hypot(br,bu),th=Math.atan2(transverse,a);
  const thp=2*Math.atan(Math.exp(-(aberr?phi:0))*Math.tan(th/2));
  const D=Math.exp(phi)*Math.cos(th/2)**2+Math.exp(-phi)*Math.sin(th/2)**2;
  const k=transverse>0?Math.sin(thp)/transverse:0;
  return {D,th,thp,f:Math.cos(thp),r:br*k,u:bu*k};
}
function boost(w,dt,X) {
  const wm=len(w);
  if(wm<1e-12)return {t:dt,X:X.slice()};
  const n=scl(w,1/wm),g=Math.cosh(wm),gb=Math.sinh(wm),par=dot(X,n);
  const tp=g*dt-gb*par,pp=g*par-gb*dt;
  return {t:tp,X:X.map((x,i)=>x+(pp-par)*n[i])};
}
function recv(p,t,along,perp,send=false) {
  const G={p,m:[1,0,0]};
  return (send?s1Ping:s1Recv)(G,t,[along,perp,0]);
}
function orbitRecv(L,a0,t,P,send=false) {
  const D=len(P),g=L.gamma;
  const lo=(t+(send?-(D+L.r):Math.max(0,D-L.r)))/g;
  const hi=(t+(send?-Math.max(0,D-L.r):D+L.r))/g;
  const first=orbitAt(L,a0,(t+(send?-D:D))/g);
  const x0=(t+(send?-1:1)*len(sub(P,first.X)))/g;
  return s1Newton(tau=>{
    const q=orbitAt(L,a0,tau),d=sub(P,q.X),R=len(d);
    return [q.t+(send?1:-1)*R-t,R>0?g*(1+(send?-1:1)*dot(q.V,d)/R):g];
  },lo,hi,x0,1e-12);
}
  function reachAS(R, n, ps) {   // end point for heading angle ps (0 = straight on) within budget n, as [along v, along e2] from you; null if n < the stopping time
    const a = R.a, c = Math.cos(ps), sn = Math.sin(ps);
    // f(t1) = t1 + phi1(t1)/a - n is non-decreasing (|d phi1/dt1| <= a): safeguarded Newton, falling back to bisection
    const fd = (t1) => { const sh = Math.sinh(a * t1), ch = Math.cosh(a * t1), ua = R.g0 * sh * c + R.gb0 * ch, us = sh * sn, um = Math.hypot(ua, us);
      const dum = um > 0 ? a * (ua * (R.g0 * ch * c + R.gb0 * sh) + us * ch * sn) / um : 0; return [t1 + Math.asinh(um) / a - n, 1 + dum / Math.sqrt(1 + um * um) / a]; };
    if (fd(0)[0] > 1e-12) return null;
    let lo = 0, hi = n, t1 = Math.min(n, Math.max(0, (n - R.phi / a) / 2));
    for (let i = 0; i < 80; i++) { const [f, d] = fd(t1); if (f < 0) lo = t1; else hi = t1; if (Math.abs(f) < 1e-12 * (1 + n) || hi - lo < 1e-13 * (1 + n)) break;
      let nt = d > 1e-9 ? t1 - f / d : NaN; if (!(nt > lo && nt < hi)) nt = (lo + hi) / 2; t1 = nt; }
    const sh = Math.sinh(a * t1), ch = Math.cosh(a * t1), k = (R.g0 - 1) * sh * c + R.gb0 * ch;
    const tM = sh / a, hM = 2 * Math.sinh(a * t1 / 2) ** 2 / a;   // (cosh - 1)/a without cancellation
    const ua = sh * c + k, us = sh * sn, um = Math.hypot(ua, us), db = um * um / (Math.sqrt(1 + um * um) + 1) / a;   // braking distance (cosh phi1 - 1)/a
    return [R.g0 * hM * c + R.gb0 * tM + (um > 0 ? ua * db / um : 0), hM * sn + (um > 0 ? us * db / um : 0)];
  }

module.exports={reachAS,trip,at,s1TauU,s1TauV,orbit,orbitAt,optics,boost,recv,orbitRecv,G_LY};
