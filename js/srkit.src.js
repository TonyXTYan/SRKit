/* srkit.js: SRKit for JavaScript. Special relativity from closed forms, computed in WebAssembly.
   This file is the source; tools/build_js.py embeds the compiled module and writes dist/srkit.js.

   One classic script, no imports and no network: it runs in a browser page, a worker, Node (require) and a
   Node vm context. The module is instantiated synchronously from the bytes embedded below, so `SRKit` is ready
   as soon as the script has run. To use a module fetched separately, call SRKit.create(bytes).

   UNITS. Natural units, c = 1: choose a unit of time and measure lengths in the distance light covers in it
   (Julian years and light-years in the game). Velocities are fractions of c; accelerations are rates of change
   of rapidity. SRKit.units holds the exact SI constants.

   SEGMENTS are plain objects, so they can be stored, copied and serialised:
     { kind: 'rest', t0, t1, P }                         at rest at P
     { kind: 'lin',  t0, t1, P0, V }                     constant velocity
     { kind: 'hyp',  t0, t1, P0, n, al, ph0 }            proper acceleration al (signed) along unit n, rapidity ph0 at t0
     { kind: 'circ', t0, t1, C, r, w, e1, e2, a0 }       circle about C in the plane (e1, e2), angle a0 + w (t - t0)
     { kind: 'null', t0, t1, P0, n, P1 }                 a light-like jump from P0 to P1
     { kind: 'steer', t0, t1, P0, U0, thrust, al }       proper acceleration al at any angle: thrust is fixed aboard, U0 = gamma v
     { kind: 'burn', t0, t1, P0, n, ph0, al0, w }        constant thrust as the mass falls; exhaust speed w, starting acceleration al0
   Positions and directions are arrays [x, y, z].

   Code AI-generated (Claude). */
(function (root) {
  'use strict';
  const EMBEDDED = '/*WASM_BASE64*/';
  const META = { version: '/*VERSION*/', sha256: '/*WASM_SHA256*/', bytes: /*WASM_BYTES*/0 };
  const ABI = 3, SEG = 21, HEADER = 24, PIECE = 24;
  const KIND = { rest: 0, lin: 1, hyp: 2, circ: 3, null: 4, steer: 5, burn: 6 };
  const KIND_NAME = ['rest', 'lin', 'hyp', 'circ', 'null', 'steer', 'burn'];
  const PHASE_NAME = ['rest', 'boost', 'coast', 'brake', 'orbit', 'jump'];
  const ZERO = [0, 0, 0];

  // base64 without atob or Buffer (neither exists in a bare vm context)
  function decode(b64) {
    const A = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/', lut = new Uint8Array(128);
    for (let i = 0; i < 64; i++) lut[A.charCodeAt(i)] = i;
    let n = b64.length; while (n && b64[n - 1] === '=') n--;
    const out = new Uint8Array(Math.floor(n * 3 / 4));
    let acc = 0, bits = 0, k = 0;
    for (let i = 0; i < n; i++) { acc = (acc << 6) | lut[b64.charCodeAt(i)]; bits += 6; if (bits >= 8) { bits -= 8; out[k++] = (acc >> bits) & 255; } }
    return out;
  }

  function create(bytes) {
    if (typeof WebAssembly === 'undefined') throw new Error('SRKit needs WebAssembly, which this environment does not provide');
    const x = new WebAssembly.Instance(new WebAssembly.Module(bytes), {}).exports, mem = x.memory;
    if (x.srkit_abi_version() !== ABI) throw new Error('SRKit: the WebAssembly module (ABI ' + x.srkit_abi_version() + ') does not match this loader (ABI ' + ABI + ')');
    let I, O, TX, calls = 0;
    function views() { const b = mem.buffer; I = new Float64Array(b, x.srkit_in(), 64); O = new Float64Array(b, x.srkit_out(), 192); TX = [0, 1, 2].map(s => new Uint8Array(b, x.srkit_text(s), 256)); }
    views();
    // when the module's memory grows, typed arrays over the old buffer are detached (their length reads 0)
    const live = () => { if (I.length === 0) views(); };
    const status = (code, what) => { if (code === 2) throw new RangeError('SRKit: invalid input to ' + what); if (code === 3) throw new RangeError('SRKit: the result of ' + what + ' does not fit a double'); return code === 0; };
    const nn = (v) => v !== v ? null : v;   // NaN means "no answer"

    // ---- segments ----
    const mk = {
      rest: (t0, t1, P) => ({ kind: 'rest', t0, t1, P }),
      lin: (t0, t1, P0, V) => ({ kind: 'lin', t0, t1, P0, V, ig: Math.sqrt(Math.max(0, 1 - (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]))) }),
      hyp: (t0, t1, P0, n, al, ph0) => ({ kind: 'hyp', t0, t1, P0, n, al, ph0 }),
      circ: (t0, t1, C, r, w, e1, e2, a0) => ({ kind: 'circ', t0, t1, C, r, w, e1, e2, a0, ig: Math.sqrt(Math.max(0, 1 - r * r * w * w)) }),
      null: (t0, t1, P0, n) => ({ kind: 'null', t0, t1, P0, n, P1: [P0[0] + n[0] * (t1 - t0), P0[1] + n[1] * (t1 - t0), P0[2] + n[2] * (t1 - t0)], ig: 0 }),
      steer: (t0, t1, P0, U0, thrust, al) => ({ kind: 'steer', t0, t1, P0, U0, thrust, al }),
      burn: (t0, t1, P0, n, ph0, al0, w) => ({ kind: 'burn', t0, t1, P0, n, ph0, al0, w }),
    };
    // write segment g into the input buffer at offset o
    function put(g, o) {
      let a, b;
      I[o + 1] = g.t0; I[o + 2] = g.t1;
      switch (g.kind) {
        case 'rest': I[o] = 0; a = g.P; I[o + 3] = a[0]; I[o + 4] = a[1]; I[o + 5] = a[2]; return;
        case 'lin': I[o] = 1; a = g.P0; b = g.V; break;
        case 'hyp': I[o] = 2; a = g.P0; b = g.n; I[o + 9] = g.al; I[o + 10] = g.ph0; break;
        case 'circ': { I[o] = 3; a = g.C; I[o + 3] = a[0]; I[o + 4] = a[1]; I[o + 5] = a[2]; I[o + 6] = g.r; I[o + 7] = g.w;
          a = g.e1; I[o + 8] = a[0]; I[o + 9] = a[1]; I[o + 10] = a[2]; a = g.e2; I[o + 11] = a[0]; I[o + 12] = a[1]; I[o + 13] = a[2]; I[o + 14] = g.a0; return; }
        case 'null': { I[o] = 4; a = g.P0; b = g.n; const p = g.P1; I[o + 9] = p[0]; I[o + 10] = p[1]; I[o + 11] = p[2]; break; }
        case 'steer': { I[o] = 5; a = g.P0; b = g.U0; const p = g.thrust; I[o + 9] = p[0]; I[o + 10] = p[1]; I[o + 11] = p[2]; I[o + 12] = g.al; break; }
        case 'burn': I[o] = 6; a = g.P0; b = g.n; I[o + 9] = g.ph0; I[o + 10] = g.al0; I[o + 11] = g.w; break;
        default: throw new Error('segment kind ' + g.kind);
      }
      I[o + 3] = a[0]; I[o + 4] = a[1]; I[o + 5] = a[2]; I[o + 6] = b[0]; I[o + 7] = b[1]; I[o + 8] = b[2];
    }
    // read a segment from the output buffer at offset o
    function take(o) {
      const k = KIND_NAME[O[o]], t0 = O[o + 1], t1 = O[o + 2], v = (i) => [O[o + i], O[o + i + 1], O[o + i + 2]];
      switch (k) {
        case 'rest': return mk.rest(t0, t1, v(3));
        case 'lin': return mk.lin(t0, t1, v(3), v(6));
        case 'hyp': return mk.hyp(t0, t1, v(3), v(6), O[o + 9], O[o + 10]);
        case 'circ': return mk.circ(t0, t1, v(3), O[o + 6], O[o + 7], v(8), v(11), O[o + 14]);
        case 'null': { const g = mk.null(t0, t1, v(3), v(6)); g.P1 = v(9); return g; }
        case 'steer': return mk.steer(t0, t1, v(3), v(6), v(9), O[o + 12]);
        case 'burn': return mk.burn(t0, t1, v(3), v(6), O[o + 9], O[o + 10], O[o + 11]);
      }
      throw new Error('segment kind code ' + O[o]);
    }
    function segAt(g, t) {
      if (g.kind === 'rest') return { X: g.P, V: ZERO };
      live(); put(g, 0); calls++; x.seg_at(t);
      return { X: [O[0], O[1], O[2]], V: [O[3], O[4], O[5]] };
    }
    function segDTau(g, t) { if (g.kind === 'rest') return t - g.t0; live(); put(g, 0); calls++; return x.seg_dtau(t); }
    function segTOf(g, s) { if (g.kind === 'rest') return g.t0 + s; live(); put(g, 0); calls++; return x.seg_t_of_tau(s); }
    function tauState() { const r = { t: O[0], X: [O[1], O[2], O[3]], V: [O[4], O[5], O[6]], phi: O[7], acc: O[8] }; if (O[9]) r.jump = true; return r; }
    function segAtTau(g, s) {
      if (g.kind === 'rest') return { t: g.t0 + s, X: g.P, V: ZERO, phi: 0, acc: 0 };
      live(); put(g, 0); calls++; x.seg_at_tau(s);
      const r = tauState(); if (g.kind === 'null') { r.X = g.P1; r.V = ZERO; } return r;
    }
    function segVmax(g) { if (g.kind === 'rest') return 0; live(); put(g, 0); calls++; return x.seg_vmax(); }
    function segRecv(g, tE, XE, lo, hi) { live(); put(g, 0); calls++; return nn(x.seg_recv(tE, XE[0], XE[1], XE[2], lo, hi)); }
    function segSend(g, tE, XE, lo, hi) { live(); put(g, 0); calls++; return nn(x.seg_send(tE, XE[0], XE[1], XE[2], lo, hi)); }
    const seg = { mk, at: segAt, dtau: segDTau, tOf: segTOf, atTau: segAtTau, vmax: segVmax, recv: segRecv, send: segSend,
      horizonMargin: (g, tE, XE, past) => { live(); put(g, 0); calls++; return nn(x.seg_horizon_margin(tE, XE[0], XE[1], XE[2], past ? 1 : 0)); },
      properDuration: (g) => { live(); put(g, 0); calls++; return x.seg_proper_duration(); } };

    // ---- worldlines owned by the module ----
    const registry = typeof FinalizationRegistry === 'function' ? new FinalizationRegistry((h) => { x.wl_free(h); }) : null;
    function build(segs, anchored, anchorT, anchorTau) {
      live(); x.wl_begin();
      for (let i = 0; i < segs.length; i++) { const g = segs[i]; put(g, 0); I[20] = g.tau0 || 0; calls++; if (!x.wl_push()) throw new RangeError('SRKit: invalid segment ' + i + ' (' + g.kind + ')'); live(); }
      calls++; const h = x.wl_finish(anchored ? 0 : 1, anchorT, anchorTau); live();
      if (!h) throw new RangeError('SRKit: a worldline needs at least one segment, in order of t0');
      return h;
    }
    class Worldline {
      /* Worldline.from(segs, anchorT = segs[0].t0, anchorTau = 0): proper time is anchorTau at anchorT.
         Call free() when done, or let the garbage collector release it. */
      constructor(h, n) { this.h = h; this.length = n; if (registry) registry.register(this, h, this); }
      static from(segs, anchorT, anchorTau = 0) { return new Worldline(build(segs, false, anchorT == null ? segs[0].t0 : anchorT, anchorTau), segs.length); }
      static fromAnchored(segs) { return new Worldline(build(segs, true, 0, 0), segs.length); }
      free() { if (this.h) { x.wl_free(this.h); if (registry) registry.unregister(this); this.h = 0; } }
      at(t) { calls++; x.wl_at(this.h, t); return { X: [O[0], O[1], O[2]], V: [O[3], O[4], O[5]], tau: O[6] }; }
      tau(t) { calls++; return x.wl_tau(this.h, t); }
      tOfTau(s) { calls++; return x.wl_t_of_tau(this.h, s); }
      atTau(s) { calls++; x.wl_at_tau(this.h, s); return tauState(); }
      find(t) { calls++; return x.wl_find(this.h, t); }
      recv(tE, XE) { calls++; return x.wl_recv(this.h, tE, XE[0], XE[1], XE[2]); }
      send(tE, XE) { calls++; return x.wl_send(this.h, tE, XE[0], XE[1], XE[2]); }
      vmax() { calls++; return x.wl_vmax(this.h); }
      continuous(positionTol = 1e-9, velocityTol = 1e-9) { calls++; return x.wl_check(this.h, positionTol, velocityTol) === 0; }
      radar(tE, XE) { calls++; x.wl_radar(this.h, tE, XE[0], XE[1], XE[2]); return { tauSend: O[0], tauRecv: O[1] }; }
      segment(i) { calls++; if (x.wl_segment(this.h, i)) return null; const g = take(0); g.tau0 = O[20]; return g; }
      /* atBatch(ts): one call for many times. Returns a Float64Array of 7 numbers per time: x, y, z, vx, vy, vz, tau. */
      atBatch(ts) {
        const n = ts.length, p = x.srkit_batch(7 * n); live();
        const B = new Float64Array(mem.buffer, p, 7 * n); B.set(ts); calls++;
        if (x.wl_at_batch(this.h, n)) throw new RangeError('SRKit: batch failed'); return B.slice();
      }
      /* coneBatch(events, send): events is a flat array of t, x, y, z. Returns the times on this worldline:
         when each event's light arrives (send false), or whose light reaches each event (send true). */
      coneBatch(events, send = false) {
        const n = events.length >> 2, p = x.srkit_batch(4 * n); live();
        const B = new Float64Array(mem.buffer, p, 4 * n); B.set(events); calls++;
        if (x.wl_cone_batch(this.h, n, send ? 1 : 0)) throw new RangeError('SRKit: batch failed');
        const out = new Float64Array(n); for (let i = 0; i < n; i++) out[i] = B[4 * i]; return out;
      }
    }

    // ---- worldlines as plain objects { segs, t0s, tau0s }: the module holds a copy and answers from it ----
    // make() registers the copy; a worldline built by hand (just { segs }) is answered segment by segment instead.
    const handles = new WeakMap();
    const wlOpt = { handles: !!registry };
    function lastLE(arr, v) { let lo = 0, hi = arr.length - 1; if (v < arr[0]) return 0; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (arr[m] <= v) lo = m; else hi = m - 1; } return lo; }
    function handleOf(W) {
      if (!wlOpt.handles) return 0;
      const e = handles.get(W); if (!e) return 0;
      const S = W.segs; return S.length === e.n && S[e.n - 1] === e.last ? e.h : 0;   // rebuilt or extended since: fall back
    }
    function wlMake(segs, anchorT, anchorTau = 0) {
      if (!segs.length) throw new Error('empty worldline');
      if (anchorT == null) anchorT = segs[0].t0;
      const W = { segs, t0s: segs.map(g => g.t0) };
      const h = build(segs, false, anchorT, anchorTau);
      for (let i = 0; i < segs.length; i++) { calls++; x.wl_segment(h, i); segs[i].tau0 = O[20]; }
      W.tau0s = segs.map(g => g.tau0);
      calls++; W.vmax = x.wl_vmax(h);
      W.hasNull = segs.some(g => g.kind === 'null');
      if (wlOpt.handles) { handles.set(W, { h, n: segs.length, last: segs[segs.length - 1] }); registry.register(W, h); } else x.wl_free(h);
      return W;
    }
    const wlFind = (W, t) => lastLE(W.t0s || W.segs.map(g => g.t0), t);
    const wlFindTau = (W, s) => lastLE(W.tau0s, s);
    function wlAt(W, t) {
      const S = W.segs; if (S.length === 1 && S[0].kind === 'rest') return { X: S[0].P, V: ZERO };
      const h = handleOf(W); if (!h) return segAt(S[wlFind(W, t)], t);
      calls++; x.wl_at(h, t); return { X: [O[0], O[1], O[2]], V: [O[3], O[4], O[5]] };
    }
    function wlTau(W, t) { const h = handleOf(W); if (h) { calls++; return x.wl_tau(h, t); } const g = W.segs[wlFind(W, t)]; return g.tau0 + segDTau(g, t); }
    function wlAtTau(W, s) {
      const h = handleOf(W); if (!h) { const g = W.segs[wlFindTau(W, s)]; return segAtTau(g, s - g.tau0); }
      calls++; x.wl_at_tau(h, s); const r = tauState(); if (r.jump) r.V = ZERO; return r;
    }
    const hasNullSeg = (W) => W.hasNull !== undefined ? W.hasNull : W.segs.some(g => g.kind === 'null');
    function wlRecv(W, tE, XE) {
      const h = handleOf(W); if (h) { calls++; return x.wl_recv(h, tE, XE[0], XE[1], XE[2]); }
      const S = W.segs, n = S.length, strict = hasNullSeg(W);
      const f = (i) => { const p = segAt(S[i], S[i].t0).X; return S[i].t0 - tE - Math.hypot(p[0] - XE[0], p[1] - XE[1], p[2] - XE[2]); };
      let a = Math.max(0, wlFind(W, tE)), b = n - 1;
      while (a < b) { const m = (a + b + 1) >> 1; if (strict ? f(m) < 0 : f(m) <= 0) a = m; else b = m - 1; }
      for (let i = a; i < n; i++) { const g = S[i], r = segRecv(g, tE, XE, i === 0 ? -Infinity : g.t0, i + 1 < n ? S[i + 1].t0 : Infinity); if (r != null) return r; }
      return Infinity;
    }
    function wlSend(W, tE, XE) {
      const h = handleOf(W); if (h) { calls++; return x.wl_send(h, tE, XE[0], XE[1], XE[2]); }
      const S = W.segs, n = S.length;
      const f = (i) => { const p = segAt(S[i], S[i].t0).X; return S[i].t0 + Math.hypot(XE[0] - p[0], XE[1] - p[1], XE[2] - p[2]) - tE; };
      let a = 0, b = wlFind(W, tE);
      while (a < b) { const m = (a + b + 1) >> 1; if (f(m) <= 0) a = m; else b = m - 1; }
      for (let i = a; i >= 0; i--) { const g = S[i], r = segSend(g, tE, XE, i === 0 ? -Infinity : g.t0, i + 1 < n ? S[i + 1].t0 : Infinity); if (r != null) return r; }
      return -Infinity;
    }
    const wl = { make: wlMake, find: wlFind, findTau: wlFindTau, at: wlAt, tau: wlTau, atTau: wlAtTau, recv: wlRecv, send: wlSend, options: wlOpt, live: () => x.wl_live() };

    // ---- itineraries ----
    const put3 = (o, v) => { I[o] = v[0]; I[o + 1] = v[1]; I[o + 2] = v[2]; };
    function takePiece(o) { const g = take(o + 3); g.tau0 = O[o + 1]; return { g, dtau: O[o + 2], phase: PHASE_NAME[O[o]] }; }
    const takeEnd = () => ({ t: O[1], tau: O[2], X: [O[3], O[4], O[5]], V: [O[6], O[7], O[8]], phi: O[9], acc: O[10] });
    function takePieces() { const n = O[0], ps = []; for (let k = 0; k < n; k++) ps.push(takePiece(HEADER + k * PIECE)); return ps; }
    function onePiece(kind, start, p1, p2, p3, p4) { calls++; x.plan_piece(kind, start.t, start.tau, p1, p2, p3, p4); return takePiece(HEADER); }
    const PHASE = { boost: 1, coast: 2, brake: 3 };
    const plan = {
      /* a start is { t, tau, X }; a piece is { g (a segment with tau0), dtau, phase } */
      rest: (start, dtau) => { live(); put3(0, start.X); return onePiece(0, start, dtau, 0, 0, 0); },
      hypPiece: (start, n, al, ph0, ph1, phase) => { live(); put3(0, start.X); put3(3, n); const p = onePiece(1, start, al, ph0, ph1, PHASE[phase] || 2); p.phase = phase; return p; },
      linPiece: (start, V, dtau) => { live(); put3(0, start.X); put3(3, V); return onePiece(2, start, dtau, 0, 0, 0); },
      nullPiece: (start, X1) => { live(); put3(0, start.X); put3(3, X1); return onePiece(3, start, 0, 0, 0, 0); },
      circPiece: (start, orb, dtau) => { live(); put3(0, start.X); put3(3, orb.C); put3(6, orb.e1); put3(9, orb.e2); return onePiece(4, start, orb.r, orb.w, orb.a0, dtau); },
      pieceEnd: (p) => { const e = segAtTau(p.g, p.dtau); e.tau = p.g.tau0 + p.dtau; return e; },
      /* line(start, X1, a, phCap, ph0, ph1, aBrake, nDefault) -> { pieces, end, n, dphi, L, phPeak, coastL } | null */
      line(start, X1, a, phCap = Infinity, ph0 = 0, ph1 = 0, aBrake = a, nDefault = null) {
        live(); put3(0, start.X); put3(3, X1); if (nDefault) put3(6, nDefault); calls++;
        if (!status(x.plan_line(start.t, start.tau, a, phCap, ph0, ph1, aBrake, nDefault ? 1 : 0), 'plan.line')) return null;
        return { pieces: takePieces(), end: takeEnd(), n: [O[11], O[12], O[13]], dphi: O[14], L: O[15], phPeak: O[16], coastL: O[17] };
      },
      brake(start, n, ph0, a) { live(); put3(0, start.X); put3(3, n); calls++; x.plan_brake(start.t, start.tau, ph0, a); return { pieces: takePieces(), end: takeEnd(), dphi: O[14] }; },
      orbitExit(start, V, a) { live(); put3(0, start.X); put3(3, V); calls++; x.plan_orbit_exit(start.t, start.tau, a); return { pieces: takePieces(), end: takeEnd(), dphi: O[14] }; },
      planeThrough(C, X, e1Default) { live(); put3(0, C); put3(3, X); if (e1Default) put3(6, e1Default); calls++; x.plan_plane_through(e1Default ? 1 : 0); return { e1: [O[0], O[1], O[2]], e2: [O[3], O[4], O[5]], nrm: [O[6], O[7], O[8]] }; },
      /* orbitEntry(start, C, r, beta, a, phCap, { e1, a0, aTangent }) -> { pieces, end, orb: { C, r, w, e1, e2, a0, beta }, dphi, fromInside } | null */
      orbitEntry(start, C, r, beta, a, phCap = Infinity, opts = {}) {
        live(); put3(0, start.X); put3(3, C); if (opts.e1) put3(6, opts.e1); calls++;
        const flags = (opts.e1 ? 1 : 0) | (opts.a0 != null ? 2 : 0) | (opts.aTangent ? 4 : 0);
        if (!status(x.plan_orbit_entry(start.t, start.tau, r, beta, a, phCap, opts.a0 || 0, opts.aTangent || 0, flags), 'plan.orbitEntry')) return null;
        return { pieces: takePieces(), end: takeEnd(), dphi: O[11], orb: { C: C.slice(), r: O[20], w: O[19], e1: [O[12], O[13], O[14]], e2: [O[15], O[16], O[17]], a0: O[18], beta: O[21] }, fromInside: !!O[22] };
      },
    };

    // ---- rides in null coordinates ----
    const tripKeys = new WeakMap();
    const keyOf = (p) => { const k = tripKeys.get(p); if (!k) throw new TypeError('SRKit: expected a ride made by SRKit.trip.plan'); return k; };
    const tripState = () => ({ tau: O[0], u: O[1], v: O[2], t: O[3], x: O[4], phi: O[5], dir: O[6], gamma: O[7], Dh: O[8], Dd: O[9] });
    const trip = {
      /* plan(d, a, cap = Infinity, stops = 0, rt = false): a ride of rest-to-rest legs at proper acceleration a */
      plan(d, a, cap = Infinity, stops = 0, rt = false) {
        const k = [d, a, cap, stops, rt ? 1 : 0]; calls++;
        status(x.trip_plan(k[0], k[1], k[2], k[3], k[4]), 'trip.plan');
        const p = { d: O[0], a: O[1], n: O[2], nOut: O[3], rt: !!O[4], legTau: O[5], tau: O[6], du: O[7], dv: O[8], uEnd: O[9], vEnd: O[10], tEnd: O[11], tOut: O[12], xEnd: O[13],
          lag: O[14], gPeak: O[15], phc: O[16], leg: { a: O[17], phc: O[18], ta: O[19], tc: O[20] } };
        tripKeys.set(p, k); return p;
      },
      at(p, tau) { const k = keyOf(p); calls++; x.trip_state(k[0], k[1], k[2], k[3], k[4], tau, 0); return tripState(); },
      atClamped(p, tau) { const k = keyOf(p); calls++; x.trip_state(k[0], k[1], k[2], k[3], k[4], tau, 1); return tripState(); },
      tauU(p, u) { const k = keyOf(p); calls++; return x.trip_tau_of(k[0], k[1], k[2], k[3], k[4], u, 0); },
      tauV(p, v) { const k = keyOf(p); calls++; return x.trip_tau_of(k[0], k[1], k[2], k[3], k[4], v, 1); },
      /* cone(p, t, along, across, send): the ship's proper time when it meets a light cone of the event */
      cone(p, t, along, across, send = false) { const k = keyOf(p); calls++; return x.trip_light_cone(k[0], k[1], k[2], k[3], k[4], t, along, across, send ? 1 : 0); },
      flipBurn(d, a) { calls++; status(x.flip_burn(d, a), 'trip.flipBurn'); return { d, a, t: O[0], tau: O[1], phiPeak: O[2], gammaPeak: O[3], lag: O[4] }; },
      /* reach(phi, a, budget, heading) -> [along, across] | null */
      reach(phi, a, budget, heading) { calls++; return status(x.reach_point(phi, a, budget, heading), 'trip.reach') ? [O[0], O[1]] : null; },
    };

    // ---- orbits, optics, rockets ----
    const takeOrbit = () => ({ r: O[0], beta: O[1], gamma: O[2], phi: O[3], omega: O[4], period: O[5], lapTau: O[6], acc: O[7], sagnac: O[8], thomas: O[9] });
    const orbit = {
      ofSpeed(r, beta) { calls++; x.orbit_from_speed(r, beta); return takeOrbit(); },
      ofAcceleration(r, acc) { calls++; status(x.orbit_from_acceleration(r, acc), 'orbit.ofAcceleration'); return takeOrbit(); },
      /* at(r, acc, a0, tau): on a circle about the origin in the xy plane */
      at(r, acc, a0, tau) { calls++; x.orbit_state(r, acc, a0, tau); return { t: O[0], X: [O[1], O[2], O[3]], V: [O[4], O[5], O[6]], heading: [O[7], O[8], O[9]], al: O[10] }; },
      cone(r, acc, a0, t, X, send = false) { calls++; return x.orbit_light_cone(r, acc, a0, t, X[0], X[1], X[2], send ? 1 : 0); },
    };
    const optics = {
      /* observe(phi, a, br, bu, aberration): direction towards the source in the observer's forward, right, up axes */
      observe(phi, a, br, bu, aberration = true) { calls++; status(x.observe(phi, a, br, bu, aberration ? 1 : 0), 'optics.observe'); return { D: O[0], th: O[1], thp: O[2], f: O[3], r: O[4], u: O[5] }; },
      playback(V, n, VB) { live(); put3(0, V); put3(3, n); put3(6, VB); calls++; return x.playback(); },
      basis(V) { live(); put3(0, V); calls++; x.basis(); return { f: [O[0], O[1], O[2]], u: [O[3], O[4], O[5]], r: [O[6], O[7], O[8]] }; },
      boost(w, dt, dX) { calls++; x.boost(w[0], w[1], w[2] || 0, dt, dX[0], dX[1], dX[2] || 0); return { t: O[3], X: [O[0], O[1], O[2]] }; },
      compose(v, u) { calls++; status(x.compose_velocities(v[0], v[1], v[2], u[0], u[1], u[2]), 'optics.compose'); return [O[0], O[1], O[2]]; },
      resolution: (aperture, wavelength, distance) => { calls++; return x.rayleigh_resolution(aperture, wavelength, distance); },
      friis: (power, dT, dR, wavelength, distance) => { calls++; return x.friis_received_power(power, dT, dR, wavelength, distance); },
      bitRate: (powerW, wavelengthM, bitsPerPhoton = 1) => { calls++; return x.photon_bit_rate(powerW, wavelengthM, bitsPerPhoton); },
      luminosity: (radiusM, temperatureK) => { calls++; return x.black_body_luminosity(radiusM, temperatureK); },
    };
    const rocket = {
      massRatio: (dphi, w = 1) => { calls++; return x.rocket_mass_ratio(dphi, w); },
      phiOf: (MR, w = 1) => { calls++; return x.rocket_rapidity(MR, w); },
      /* energy(dphi, w, efficiency): propellant energy per unit of dry rest energy; times m c^2 for joules */
      energy: (dphi, w = 1, efficiency = 1) => { calls++; return x.rocket_energy(dphi, w, efficiency); },
      sailPhi: (energyPerRestEnergy) => { calls++; return x.sail_rapidity(energyPerRestEnergy); },
    };
    const kepler = (M, e) => { calls++; return x.kepler_eccentric_anomaly(M, e); };

    // ---- exact time: whole milliseconds as decimal strings or BigInt, 384 or 512 bits wide ----
    function exact(width) {
      const f = (name) => x['time' + width + '_' + name];
      const putText = (slot, v) => { const s = String(v); if (s.length > 256) throw new RangeError('SRKit: number too long'); live(); for (let i = 0; i < s.length; i++) TX[slot][i] = s.charCodeAt(i); return s.length; };
      const takeText = (len, what) => { if (len < 0) throw new RangeError('SRKit: ' + what + (len === -3 ? ' is out of range' : ' has invalid input')); live(); let s = ''; for (let i = 0; i < len; i++) s += String.fromCharCode(TX[2][i]); return s; };
      const like = (ref, s) => typeof ref === 'bigint' ? BigInt(s) : s;
      return {
        bits: width,
        /* add(t, d): timestamp plus signed duration, both in ms */
        add: (t, d) => { calls++; return like(t, takeText(f('add')(putText(0, t), putText(1, d)), 'time.add')); },
        /* diff(a, b): a - b as a signed duration in ms */
        diff: (a, b) => { calls++; return like(a, takeText(f('diff')(putText(0, a), putText(1, b)), 'time.diff')); },
        cmp: (a, b) => { calls++; const r = f('cmp')(putText(0, a), putText(1, b)); if (r === -2) throw new RangeError('SRKit: time.cmp has invalid input'); return r; },
        /* mod(t, periodMs): the exact phase of t in a period of whole milliseconds */
        mod: (t, periodMs) => { calls++; return nn(f('mod')(putText(0, t), periodMs)); },
        /* seconds(a, b): a - b in seconds, subtracted exactly first; null beyond 2^53 ms */
        seconds: (a, b) => { calls++; return nn(f('seconds')(putText(0, a), putText(1, b))); },
        /* lightTicks(dX): light-travel time over dX metres, rounded up to whole ms */
        lightTicks: (dX) => { calls++; return takeText(f('light_ticks')(dX[0], dX[1], dX[2]), 'time.lightTicks'); },
        firstArrival: (emission, dX) => { calls++; return like(emission, takeText(f('first_arrival')(putText(0, emission), dX[0], dX[1], dX[2]), 'time.firstArrival')); },
        /* arrived(emission, now, dX): has the signal crossed dX metres by tick now? exact */
        arrived: (emission, now, dX) => { calls++; const r = f('arrived')(putText(0, emission), putText(1, now), dX[0], dX[1], dX[2]); if (r < 0) throw new RangeError('SRKit: time.arrived has invalid input'); return r === 1; },
        /* classify(dt, dX): 'timelike' | 'null' | 'spacelike', exact */
        classify: (dt, dX) => { calls++; const r = f('classify')(putText(0, dt), dX[0], dX[1], dX[2]); if (r === -2) throw new RangeError('SRKit: time.classify has invalid input'); return r === 1 ? 'timelike' : r === 0 ? 'null' : 'spacelike'; },
      };
    }
    const U = (i) => x.unit_constant(i);
    const units = Object.freeze({ C: U(0), YEAR: U(1), LY: U(2), AU: U(3), PC: U(4), G0: U(5), H: U(6), HBAR: U(7), KB: U(8), G: U(9), SIGMA: U(10), L_SUN: U(11), R_SUN: U(12), GM_SUN: U(13),
      G_PER_YEAR: U(14), AU_LY: U(15), PC_LY: U(16) });

    /*EXTRA_BINDINGS*/

    return Object.freeze({
      info: Object.freeze({ backend: 'srkit-wasm', version: META.version, abi: ABI, sha256: META.sha256, bytes: META.bytes || bytes.length }),
      calls: () => calls, units, seg, wl, Worldline, plan, trip, orbit, optics, rocket, kepler,
      time: Object.freeze(Object.assign(exact(384), { u384: exact(384), u512: exact(512) })),
      /*EXTRA_EXPORTS*/
      raw: x, create,
    });
  }

  let api = null;
  if (EMBEDDED.length > 64) api = create(decode(EMBEDDED));
  else api = Object.freeze({ create });   // no module embedded: call SRKit.create(bytes)
  if (typeof module === 'object' && module && module.exports) module.exports = api;
  root.SRKit = api;
})(typeof window !== 'undefined' ? window : typeof globalThis !== 'undefined' ? globalThis : this);
