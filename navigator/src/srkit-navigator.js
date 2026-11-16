/* The Relativistic Navigator's physics, on SRKit. AI-generated code.

   The page calls SRKitNavigator; this file maps those calls onto SRKit's general API (window.SRKit, the Rust engine
   compiled to WebAssembly and embedded above this script by build_srkit.py). Nothing here computes physics: it
   converts the page's units (drives in g, light-years and Julian years with c = 1) and shapes results into the field
   names the page was written against. There is no JavaScript fallback: without the engine the page says so and stops. */
(function (root) {
  'use strict';
  const K = root.SRKit;
  if (!K || !K.trip) {
    if (typeof document !== 'undefined' && document.body) {
      const message = document.createElement('p'); message.setAttribute('role', 'alert');
      message.textContent = 'The physics engine could not start. Open this page in a browser with WebAssembly enabled.';
      document.body.prepend(message);
    }
    throw new Error('SRKit WebAssembly startup failed: the engine did not load');
  }
  const G = K.units.G_PER_YEAR;   // 1 g as a rate of change of rapidity per year: 1.0323
  const unit = (v) => { const n = Math.hypot(v[0], v[1], v[2] || 0); if (!(n > 0) || !isFinite(n)) throw new RangeError('SRKit: invalid input to travellerEmission (heading)'); return [v[0] / n, v[1] / n, (v[2] || 0) / n]; };
  const api = {
    info: K.info,
    calls: K.calls,
    /* a ride of rest-to-rest legs: distance d (ly), drive g, peak gamma capped at cap, evenly spaced stops, and back if rt */
    trip: (d, g, cap = Infinity, stops = 0, rt = false) => K.trip.plan(d, g * G, cap, stops, rt),
    /* the ship at ship time tau; it waits at rest before launch and after arrival */
    at: (p, tau) => K.trip.at(p, tau),
    tauU: (p, u) => K.trip.tauU(p, u),
    tauV: (p, v) => K.trip.tauV(p, v),
    /* a circle of radius r (ly) flown at a drive of g */
    orbit(r, g) {
      const o = K.orbit.ofAcceleration(r, g * G);
      return { r, a: o.acc, gamma: o.gamma, beta: o.beta, phi: o.phi, lapT: o.period, lapTau: o.lapTau, omega: o.omega, sagnac: o.sagnac };
    },
    orbitAt: (p, a0, tau) => K.orbit.at(p.r, p.a, a0, tau),
    observe: (phi, a, br, bu, aberr = true) => K.optics.observe(phi, a, br, bu, aberr),
    boost: (w, dt, dX) => K.optics.boost(w, dt, dX),
    /* ship time at which the ride meets a light cone of the event (t, along, across): received, or sent if send */
    recv: (p, t, a, b, send = false) => K.trip.cone(p, t, a, b, send),
    orbitRecv: (p, a0, t, P, send = false) => K.orbit.cone(p.r, p.a, a0, t, P, send),
    reach: (phi, a, budget, heading) => K.trip.reach(phi, a, budget, heading),
    /* when a test traveller emitted the light now seen: it waits at the origin until t = 0, then coasts at beta along
       heading. The observer is at relativePosition (from the traveller's launch point) at time t. A worldline of two
       closed-form segments, and one light-cone query. */
    travellerEmission(beta, heading, t, relativePosition) {
      if (!(beta >= 0 && beta < 1) || !isFinite(t) || !relativePosition.every(isFinite)) throw new RangeError('SRKit: invalid input to travellerEmission');
      const n = unit(heading);
      const W = { segs: [K.seg.mk.rest(-Infinity, 0, [0, 0, 0]), K.seg.mk.lin(0, Infinity, [0, 0, 0], [n[0] * beta, n[1] * beta, n[2] * beta])], t0s: [-Infinity, 0] };
      return K.wl.send(W, t, [relativePosition[0], relativePosition[1], relativePosition[2] || 0]);
    },
  };
  root.SRKitNavigator = Object.freeze(api);
})(typeof window !== 'undefined' ? window : globalThis);
