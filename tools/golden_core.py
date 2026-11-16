#!/usr/bin/env python3
"""Golden vectors for srkit-core, from first principles at 60 significant digits.

Each row is `case;inputs;outputs`: inputs are binary64 values written so that they read back exactly, outputs are
the exact mathematical result of those inputs to 25 digits. The Rust test `tests/golden.rs` evaluates the same
cases and compares to a stated relative tolerance. Nothing here shares code with the Rust crate: every case is its
textbook formula (or a quadrature, or a root found by mpmath) and is independent of how the crate computes it.

Run from the repository root:
    python srkit/tools/golden_core.py > srkit/crates/srkit-core/tests/golden/core.txt
The output is deterministic (fixed seed). Changing a golden value is a breaking change: see the README.
"""
import random
import sys

import mpmath as mp

mp.mp.dps = 60
rng = random.Random(20261116)
out = []


def f64(x):
    return float(x)


def emit(case, inputs, outputs):
    out.append(case + ';' + ','.join(repr(f64(v)) for v in inputs) + ';' + ','.join(mp.nstr(mp.mpf(v), 25) for v in outputs))


def logu(lo, hi):
    return f64(10 ** rng.uniform(lo, hi))


# --- hyperbolic segment, by coordinate time: rapidity phi0 along n, proper acceleration alpha, after dt ---
def hyp_at(phi0, alpha, dt):
    p, a, d = mp.mpf(phi0), mp.mpf(alpha), mp.mpf(dt)
    sh = mp.sinh(p) + a * d
    phi = mp.asinh(sh)
    return [(mp.cosh(phi) - mp.cosh(p)) / a, mp.tanh(phi), (phi - p) / a]


for _ in range(60):   # boosting from rest or from speed, any scale from light-seconds to the far side of the universe
    phi0, alpha, dt = rng.choice([0.0, 0.0, logu(-6, 1.2)]), logu(-3, 2), logu(-9, 12)
    emit('hyp_at', [phi0, alpha, dt], hyp_at(phi0, alpha, dt))
for _ in range(30):   # braking, kept clear of the turnaround (where coordinate time cannot resolve the state)
    phi0, alpha = logu(-2, 1.2), logu(-3, 2)
    frac = rng.uniform(0.05, 0.9)
    dt = f64(mp.sinh(phi0) * frac / alpha)
    emit('hyp_at', [phi0, -alpha, dt], hyp_at(phi0, -alpha, dt))

# --- hyperbolic segment, by proper time ---
for _ in range(60):
    phi0, alpha, s = rng.choice([0.0, logu(-6, 1.2)]) * rng.choice([1, 1, -1]), logu(-3, 2) * rng.choice([1, -1]), logu(-9, 1.5)
    p, a, t = mp.mpf(phi0), mp.mpf(alpha), mp.mpf(s)
    phi = p + a * t
    emit('hyp_tau', [phi0, alpha, s], [(mp.sinh(phi) - mp.sinh(p)) / a, (mp.cosh(phi) - mp.cosh(p)) / a, phi])

# --- flip-and-burn from rest to rest ---
for _ in range(60):
    d, a = logu(-9, 11), logu(-4, 3)
    D, A = mp.mpf(d), mp.mpf(a)
    phi = mp.acosh(1 + A * D / 2)
    t = mp.sqrt(D * D + 4 * D / A)
    emit('flipburn', [d, a], [t, 2 * phi / A, phi, 1 + A * D / 2, t - D])


# --- a ride in null coordinates: legs with a cap, stops and a return ---
def trip_ref(d, a, cap, stops, rt):
    D, A = mp.mpf(d), mp.mpf(a)
    n_out = stops + 1
    n = 2 * n_out if rt else n_out
    L = D / n_out
    gpk = 1 + A * L / 2
    if mp.mpf(cap) < gpk:
        phc = mp.acosh(mp.mpf(cap))
        tc = (L - 2 * (mp.mpf(cap) - 1) / A) / mp.sinh(phc)
    else:
        phc, tc = mp.acosh(gpk), mp.mpf(0)
    ta = phc / A
    du = 2 * (1 - mp.exp(-phc)) / A + mp.exp(-phc) * tc
    dv = 2 * (mp.exp(phc) - 1) / A + mp.exp(phc) * tc
    return dict(A=A, n=n, n_out=n_out, phc=phc, ta=ta, tc=tc, du=du, dv=dv, leg=2 * ta + tc, rt=rt)


def trip_state(T, tau):
    A, phc, ta, tc = T['A'], T['phc'], T['ta'], T['tc']
    i = min(T['n'] - 1, int(mp.floor(tau / T['leg'])))
    s = tau - i * T['leg']
    if s <= ta:
        phi = A * s
        lu, lv = (1 - mp.exp(-phi)) / A, (mp.exp(phi) - 1) / A
    elif s <= ta + tc:
        phi = phc
        lu = (1 - mp.exp(-phc)) / A + mp.exp(-phc) * (s - ta)
        lv = (mp.exp(phc) - 1) / A + mp.exp(phc) * (s - ta)
    else:
        phi = phc - A * (s - ta - tc)
        lu = (1 - mp.exp(-phc)) / A + mp.exp(-phc) * tc + (mp.exp(-phi) - mp.exp(-phc)) / A
        lv = (mp.exp(phc) - 1) / A + mp.exp(phc) * tc + (mp.exp(phc) - mp.exp(phi)) / A
    back = T['rt'] and i >= T['n_out']
    if back:
        j = i - T['n_out']
        return T['n_out'] * T['du'] + j * T['dv'] + lv, T['n_out'] * T['dv'] + j * T['du'] + lu, phi
    return i * T['du'] + lu, i * T['dv'] + lv, phi


for _ in range(60):
    d, a = logu(-2, 6.5), logu(-2, 2)
    cap = rng.choice([float('inf'), float('inf'), logu(0.05, 4)])
    stops, rt = rng.choice([0, 0, 1, 3]), rng.choice([0, 0, 1])
    T = trip_ref(d, a, cap, stops, rt)
    n = T['n']
    emit('trip', [d, a, cap, stops, rt], [n * T['leg'], n * (T['du'] + T['dv']) / 2, T['du'], T['dv'], mp.cosh(T['phc']), T['n_out'] * T['du']])
    # a state strictly inside a leg (not within rounding of a leg boundary, where the leg index is a convention)
    frac = (rng.randrange(n) + rng.uniform(0.02, 0.98)) / n
    tau = f64(n * T['leg'] * frac)
    u, v, phi = trip_state(T, mp.mpf(tau))
    emit('trip_state', [d, a, cap, stops, rt, tau], [u, v, phi])

# --- aberration and Doppler: observer at rapidity phi, source at angle theta from the direction of motion ---
for _ in range(60):
    phi, theta = logu(-8, 1.3) * rng.choice([1, 1, 1, -1]), rng.uniform(1e-6, 3.14159)
    p, th = mp.mpf(phi), mp.mpf(theta)
    D = mp.cosh(p) + mp.sinh(p) * mp.cos(th)
    emit('observe', [phi, theta], [D, 2 * mp.atan(mp.exp(-p) * mp.tan(th / 2))])


# --- boost by a rapidity vector ---
def vec(scale):
    return [rng.uniform(-1, 1) * scale for _ in range(3)]


for _ in range(40):
    w, dt, dx = vec(logu(-6, 1)), rng.uniform(-1, 1) * logu(-3, 6), vec(logu(-3, 6))
    W = mp.matrix(w)
    phi = mp.norm(W)
    n = W / phi
    X = mp.matrix(dx)
    par = (n.T * X)[0]
    g, gb = mp.cosh(phi), mp.sinh(phi)
    tp = g * mp.mpf(dt) - gb * par
    pp = g * par - gb * mp.mpf(dt)
    Xp = X + (pp - par) * n
    emit('boost', w + [dt] + dx, [tp, Xp[0], Xp[1], Xp[2]])

# --- thrust at an angle fixed aboard: U(q) = U0 cosh q + N sinh q ---
for _ in range(60):
    speed = rng.choice([0.0, rng.uniform(0.01, 0.99), 1 - logu(-6, -2)])
    v0 = [0.0, 0.0, 0.0]
    if speed:
        d = mp.matrix(vec(1)); d = d / mp.norm(d)
        v0 = [f64(d[i] * speed) for i in range(3)]
    thrust, alpha, dt = vec(1), logu(-2, 1.5), logu(-4, 3)
    V0 = mp.matrix(v0)
    g0 = 1 / mp.sqrt(1 - (V0.T * V0)[0])
    U0 = g0 * V0
    nn = mp.matrix(thrust); nn = nn / mp.norm(nn)
    nt = (U0.T * nn)[0]
    Ns = nn + U0 * (nt / (g0 + 1))
    A = mp.mpf(alpha)
    T = A * mp.mpf(dt)
    # g0 sinh q + nt (cosh q - 1) = T, as a quadratic in e^q
    E = ((T + nt) + mp.sqrt((T + nt) ** 2 + g0 * g0 - nt * nt)) / (g0 + nt)
    q = mp.log(E)
    X = (U0 * mp.sinh(q) + Ns * (mp.cosh(q) - 1)) / A
    U = U0 * mp.cosh(q) + Ns * mp.sinh(q)
    gam = g0 * mp.cosh(q) + nt * mp.sinh(q)
    emit('steered', v0 + thrust + [alpha, dt], [X[0], X[1], X[2], U[0] / gam, U[1] / gam, U[2] / gam, q / A])

# --- constant thrust while the mass falls: quadrature of dt = cosh(phi) ds, dx = sinh(phi) ds ---
for _ in range(40):
    phi0, w = rng.choice([0.0, logu(-3, 0.5)]), rng.choice([1.0, rng.uniform(0.01, 0.99)])
    alpha0 = logu(-2, 1) * rng.choice([1, 1, -1])
    tb = mp.mpf(w) / abs(mp.mpf(alpha0))
    s = f64(tb * rng.uniform(0.001, 0.95))
    sg = 1 if alpha0 > 0 else -1
    phi = lambda x: mp.mpf(phi0) - sg * mp.mpf(w) * mp.log(1 - x / tb)
    S = mp.mpf(s)
    emit('burn', [phi0, alpha0, w, s], [mp.quad(lambda x: mp.cosh(phi(x)), mp.linspace(0, S, 9)), mp.quad(lambda x: mp.sinh(phi(x)), mp.linspace(0, S, 9)), phi(S), 1 - S / tb])

# --- Kepler's equation ---
for _ in range(40):
    M, e = rng.uniform(-6.3, 6.3), rng.uniform(0, 0.95)
    emit('kepler', [M, e], [mp.findroot(lambda E: E - mp.mpf(e) * mp.sin(E) - mp.mpf(M), mp.mpf(M))])

# --- reach: burn at a fixed heading for t1, then brake to rest along the velocity; t1 + phi1/a = budget ---
for _ in range(40):
    phi, a = rng.choice([0.0, logu(-3, 0.8)]), logu(-1, 1)
    budget = f64(phi / a + logu(-2, 1.2))
    heading = rng.uniform(0, 3.14159)
    P, A, B, H = mp.mpf(phi), mp.mpf(a), mp.mpf(budget), mp.mpf(heading)
    g0, gb0, c, sn = mp.cosh(P), mp.sinh(P), mp.cos(H), mp.sin(H)

    def after(t1):
        sh, ch = mp.sinh(A * t1), mp.cosh(A * t1)
        return g0 * sh * c + gb0 * ch, sh * sn   # proper velocity along and across after the burn

    t1 = mp.findroot(lambda t: t + mp.asinh(mp.hypot(*after(t))) / A - B, (0, B), solver='anderson', tol=1e-45, maxsteps=200)
    sh, ch = mp.sinh(A * t1), mp.cosh(A * t1)
    ua, us = after(t1)
    um = mp.hypot(ua, us)
    db = (mp.sqrt(1 + um * um) - 1) / A
    along = g0 * (ch - 1) / A * c + gb0 * sh / A + (ua * db / um if um else 0)
    across = (ch - 1) / A * sn + (us * db / um if um else 0)
    emit('reach', [phi, a, budget, heading], [along, across])

sys.stdout.write('\n'.join(out) + '\n')
