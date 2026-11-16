"""Independent high-precision reference values for srkit-gr.

Writes tests/data/crosscheck.csv, which tests/crosscheck.rs reads. Each row names a public function, its
arguments (exact binary64 values, printed so that they parse back to the same bits) and the value mpmath gives
at 120 significant digits for those exact arguments and the exact binary64 constants of srkit_core::units.

The references do not copy the crate's formulas where another route exists:

  * Kerr clock rates, frame dragging and the energy, angular momentum and dt/dtau of circular orbits come from
    the Boyer-Lindquist metric components (g_tt, g_tphi, g_phiphi), not from the rearranged closed forms the
    crate evaluates. The horizon angular velocity and the irreducible mass come from the horizon's
    g_phiphi and area.
  * The characteristic radii are tied to their defining conditions: the photon orbit and the marginally bound
    orbit are roots found numerically (of the vanishing denominator, and of E = 1); the ISCO is the
    Bardeen-Press-Teukolsky formula in its original form, checked against the condition
    x^2 - 6x + 8a sqrt(x) - 3a^2 = 0 to 50 digits. Horizons and the static limit are the roots of their
    quadratics, which is also what the crate computes: those rows test the arithmetic only.
  * Shapiro delays are numerical integrals of 2GM/(c^3 r) along the straight path.
  * The deflection angle is the integral of the transverse acceleration, the Einstein angle a root of the lens
    equation, the binding energy an integral over shells, the Hawking photon-black-body power sigma A T^4, the
    lifetime an integral of dM/(dM/dt), and the Hohmann burns come from the vis-viva equation.
  * Where a formula is a single product or quotient (a radius k GM/c^2, the Hawking temperature, a circular
    speed, a clock-rate offset) there is no second route, and the row tests the arithmetic only.

Run (from the repository root; needs mpmath only):

    python srkit/crates/srkit-gr/tools/crosscheck.py            # rewrite the CSV
    python srkit/crates/srkit-gr/tools/crosscheck.py --check    # fail if the checked-in CSV is out of date

The constants below must equal srkit-core/src/units.rs. They are written into the CSV as `constant` rows and
the Rust test compares them bit for bit, so a change of a constant fails the test until this file is rerun.
"""
import math
import random
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 120

# ---- constants: must match srkit_core::units ----
C = 299_792_458.0
G = 6.674_30e-11
HBAR = 1.054_571_817_646_156_5e-34
KB = 1.380_649e-23
GM_SUN = 1.327_124_4e20
M_SUN = GM_SUN / G                      # the same binary64 division as the Rust constant
AU = 149_597_870_700.0
R_SUN = 6.957e8
JULIAN_YEAR = 31_557_600.0
CONSTANTS = [('C', C), ('G_NEWTON', G), ('HBAR', HBAR), ('BOLTZMANN', KB), ('GM_SUN', GM_SUN), ('M_SUN', M_SUN),
             ('AU', AU), ('R_SUN', R_SUN), ('JULIAN_YEAR', JULIAN_YEAR)]

c, g, hbar, kb = (mp.mpf(v) for v in (C, G, HBAR, KB))
pi = mp.pi

OUT = Path(__file__).resolve().parent.parent / 'tests' / 'data' / 'crosscheck.csv'
rows = []
random.seed(20261116)


def emit(function, args, expected, tol):
    """One row. `args` are binary64 values (or 'inf'); `expected` an mpf; `tol` the relative tolerance."""
    cells = [repr(float(a)) for a in args] + [''] * (4 - len(args))
    rows.append(','.join([function] + cells + [mp.nstr(expected, 25, min_fixed=0, max_fixed=0), repr(tol)]))


# Default tolerance: "within a few ULP". One ULP is at most 2.2e-16 relative.
ULPS = 2e-15


def grav_radius(mass):
    return g * mp.mpf(mass) / c ** 2


# ======================= top level =======================
MASSES = [1.0, 1e9, 5.9722e24, M_SUN, 10.0 * M_SUN, 4.297e6 * M_SUN, 1e11 * M_SUN]
for M in MASSES:
    emit('gravitational_parameter', [M], g * mp.mpf(M), ULPS)
    emit('gravitational_radius', [M], grav_radius(M), ULPS)
    emit('gravitational_time', [M], grav_radius(M) / c, ULPS)
    emit('mass_from_gravitational_parameter', [M], mp.mpf(M) / g, ULPS)

# ======================= Schwarzschild =======================
# radii in units of GM/c^2, converted to binary64 metres; the reference uses the rounded metres exactly
SCHW_X = [2.0 + 1e-12, 2.0 + 1e-9, 2.000001, 2.001, 2.1, 2.5, 3.0 + 1e-9, 3.001, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 30.0,
          100.0, 1e3, 1e5, 1e8, 1e12]
for M in [M_SUN, 5.9722e24]:
    m = grav_radius(M)
    emit('schwarzschild::schwarzschild_radius', [M], 2 * m, ULPS)
    emit('schwarzschild::photon_sphere_radius', [M], 3 * m, ULPS)
    emit('schwarzschild::marginally_bound_radius', [M], 4 * m, ULPS)
    emit('schwarzschild::isco_radius', [M], 6 * m, ULPS)
    radii = [float(m * x) for x in SCHW_X]
    for r in radii:
        rr = mp.mpf(r)
        lapse2 = 1 - 2 * m / rr
        if lapse2 > 0:
            # Near the horizon the result is as well conditioned as r - r_s, which the crate keeps exact; the
            # rounding of r_s itself (a few ULP of r_s) is what limits it.
            cond = max(1.0, float(1 / lapse2))
            emit('schwarzschild::static_dtau_dt', [M, r], mp.sqrt(lapse2), ULPS * cond)
        if rr > 3 * m:
            omega2 = g * mp.mpf(M) / rr ** 3
            # u^t from the metric: -1 = g_tt (u^t)^2 c^2 + r^2 Omega^2 (u^t)^2  (equatorial circular orbit)
            ut = 1 / mp.sqrt(lapse2 - rr ** 2 * omega2 / c ** 2)
            cond = max(1.0, float(1 / (1 - 3 * m / rr)))
            emit('schwarzschild::circular_orbit_dt_dtau', [M, r], ut, ULPS * cond)
            emit('schwarzschild::circular_orbit_angular_velocity', [M, r], mp.sqrt(omega2), ULPS)
    emit('schwarzschild::static_dtau_dt', [M, math.inf], mp.mpf(1), ULPS)
    # redshift between pairs of radii, including to infinity and nearly equal radii
    pairs = [(radii[i], radii[j]) for i in range(1, len(radii), 3) for j in range(2, len(radii), 4) if i != j]
    pairs += [(radii[4], math.inf), (radii[12], math.inf), (math.inf, radii[9]), (radii[16], radii[16] * (1 + 1e-9))]
    for re, ro in pairs:
        ae = 1 - 2 * m / mp.mpf(re) if re != math.inf else mp.mpf(1)
        ao = 1 - 2 * m / mp.mpf(ro) if ro != math.inf else mp.mpf(1)
        if ae <= 0 or ao <= 0:
            continue
        z = mp.sqrt(ao / ae) - 1
        cond = max(1.0, float(1 / ae), float(1 / ao))
        emit('schwarzschild::gravitational_redshift', [M, re, ro], z, ULPS * cond)

# ======================= Kerr =======================
SPINS = [0.0, 1e-30, 1e-12, 1e-6, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.998, 0.9999, 1 - 1e-6, 1 - 1e-9,
         1 - 1e-12, 1 - 1e-14, 1.0]
SPINS += [random.random() for _ in range(3)] + [1 - 10 ** random.uniform(-14, -1) for _ in range(3)]
PRO, RETRO = 1.0, -1.0


def metric(x, a, theta):
    """Boyer-Lindquist g_tt, g_tphi, g_phiphi with M = 1 (signature -+++)."""
    s2, c2 = mp.sin(theta) ** 2, mp.cos(theta) ** 2
    sigma = x * x + a * a * c2
    gtt = -(1 - 2 * x / sigma)
    gtp = -2 * a * x * s2 / sigma
    gpp = (x * x + a * a + 2 * a * a * x * s2 / sigma) * s2
    return gtt, gtp, gpp


def circular(x, a):
    """u^t, E, L, Omega of the equatorial circular geodesic at x (M = 1), from the metric and Omega."""
    omega = 1 / (x ** mp.mpf('1.5') + a)
    gtt, gtp, gpp = metric(x, a, pi / 2)
    ut = 1 / mp.sqrt(-(gtt + 2 * gtp * omega + gpp * omega ** 2))
    return ut, -(gtt + gtp * omega) * ut, (gtp + gpp * omega) * ut, omega


def isco(a):
    """ISCO radius (M = 1) for signed spin a: the root of x^2 - 6x + 8a sqrt(x) - 3a^2 outside the photon orbit."""
    if abs(a) == 1:
        return mp.mpf(1) if a > 0 else mp.mpf(9)
    z1 = 1 + mp.cbrt(1 - a * a) * (mp.cbrt(1 + a) + mp.cbrt(1 - a))
    z2 = mp.sqrt(3 * a * a + z1 * z1)
    x = 3 + z2 - mp.sign(a) * mp.sqrt((3 - z1) * (3 + z1 + 2 * z2))
    residual = x * x - 6 * x + 8 * a * mp.sqrt(x) - 3 * a * a
    assert abs(residual) < mp.mpf('1e-50'), (a, residual)
    return x


def photon(a):
    """Photon orbit radius (M = 1): the root of x^(3/2) - 3 sqrt(x) + 2a between the horizon and 4."""
    if abs(a) == 1:
        return mp.mpf(1) if a > 0 else mp.mpf(4)
    start = 2 * (1 + mp.cos(mp.mpf(2) / 3 * mp.acos(-a)))
    x = mp.findroot(lambda t: t ** mp.mpf('1.5') - 3 * mp.sqrt(t) + 2 * a, start)
    assert abs(x - start) < mp.mpf('1e-40')
    return x


def marginally_bound(a):
    """Marginally bound radius (M = 1): the circular orbit with E = 1."""
    if a == 1:
        return mp.mpf(1)
    start = 2 - a + 2 * mp.sqrt(1 - a)
    x = mp.findroot(lambda t: circular(t, a)[1] - 1, start)
    assert abs(x - start) < mp.mpf('1e-40')
    return x


KERR_MASS = 10.0 * M_SUN
m = grav_radius(KERR_MASS)

# The hole itself: both signs of the spin for a few values, the rest positive.
for spin0 in SPINS:
    for sign in ((1.0, -1.0) if spin0 in (0.3, 0.998, 1.0) else (1.0,)):
        spin = sign * spin0
        a = mp.mpf(spin)
        root = mp.sqrt(1 - a * a)
        xp, xm = 1 + root, 1 - root                 # roots of Delta = x^2 - 2x + a^2
        emit('kerr::outer_horizon_radius', [KERR_MASS, spin], m * xp, ULPS)
        emit('kerr::inner_horizon_radius', [KERR_MASS, spin], m * xm, ULPS)
        emit('kerr::horizon_angular_velocity', [KERR_MASS, spin], a / (xp * xp + a * a) * c / m, ULPS)
        m_irr = mp.sqrt((xp * xp + a * a) / 4) * mp.mpf(KERR_MASS)   # horizon area 4 pi (r+^2 + a^2) = 16 pi M_irr^2
        emit('kerr::irreducible_mass', [KERR_MASS, spin], m_irr, ULPS)
        emit('kerr::extractable_energy_fraction', [spin], 1 - m_irr / mp.mpf(KERR_MASS), ULPS)
        for theta in (0.3, math.pi / 2):
            xe = 1 + mp.sqrt(1 - a * a * mp.cos(mp.mpf(theta)) ** 2)        # outer root of g_tt = 0
            emit('kerr::ergosphere_radius', [KERR_MASS, spin, theta], m * xe, ULPS)

# Observers: radii as multiples of the horizon (ZAMO) and of the static limit (static observer).
for spin in (0.0, 0.5, 0.9, 0.998, 1 - 1e-9, 1.0, -0.7):
    a = mp.mpf(spin)
    xp = 1 + mp.sqrt(1 - a * a)
    for theta in ((1.0,) if spin < 0 else (0.0, 1.0, math.pi / 2)):
        th = mp.mpf(theta)
        xe = 1 + mp.sqrt(1 - a * a * mp.cos(th) ** 2)
        for k in (1 + 1e-9, 1.01, 2.5, 1e3):
            r = float(m * xp * k)
            x = mp.mpf(r) / m
            gtt, gtp, gpp = metric(x, a, th)
            if theta == 0.0:
                # on the axis g_phiphi = 0: the limit of -g_tphi/g_phiphi is 2ax/(x^2 + a^2)^2
                omega = 2 * a * x / (x * x + a * a) ** 2
                lapse2 = (x * x - 2 * x + a * a) / (x * x + a * a)
            else:
                omega = -gtp / gpp
                lapse2 = -gtt + gtp * gtp / gpp
            # close to the horizon the lapse is as well conditioned as r - r+
            cond = max(1.0, float(x / (x - xp)))
            emit('kerr::zamo_angular_velocity', [KERR_MASS, spin, r, theta], omega * c / m, ULPS * 2)
            emit('kerr::zamo_dtau_dt', [KERR_MASS, spin, r, theta], mp.sqrt(lapse2), ULPS * cond)
            rs = float(m * xe * k)
            xs = mp.mpf(rs) / m
            cond = max(1.0, float(xs / (xs - xe)))
            emit('kerr::static_dtau_dt', [KERR_MASS, spin, rs, theta], mp.sqrt(-metric(xs, a, th)[0]), ULPS * cond)

# Orbits depend on |spin| and the sense.
for spin0 in SPINS:
    for sense in ((PRO,) if spin0 == 0.0 else (PRO, RETRO)):
        a = mp.mpf(spin0) * int(sense)
        x_isco, x_ph, x_mb = isco(a), photon(a), marginally_bound(a)
        emit('kerr::isco_radius', [KERR_MASS, spin0, sense], m * x_isco, ULPS)
        emit('kerr::photon_orbit_radius', [KERR_MASS, spin0, sense], m * x_ph, ULPS)
        emit('kerr::marginally_bound_radius', [KERR_MASS, spin0, sense], m * x_mb, ULPS)
        if a == 1:
            continue        # extremal prograde ISCO: dt/dtau is infinite; the limits are tested in tests/limits.rs
        ut, energy, ang, _ = circular(x_isco, a)
        # dt/dtau on the prograde ISCO grows as (1 - a)^(-1/3) near the extremal limit and is only as good as
        # x_isco - 1, which binary64 resolves to 2.2e-16/(x_isco - 1)
        cond = max(1.0, float(1 / (x_isco - 1)))
        emit('kerr::isco_specific_angular_momentum', [KERR_MASS, spin0, sense], ang * m * c, ULPS)
        emit('kerr::isco_dt_dtau', [spin0, sense], ut, ULPS * 2 * cond)
        emit('kerr::isco_specific_energy', [spin0, sense], energy, ULPS)
        emit('kerr::accretion_efficiency', [spin0, sense], 1 - energy, ULPS)

# Circular orbits at multiples of the photon orbit radius, and on the ISCO itself.
for spin0 in (0.0, 0.3, 0.9, 0.998, 1 - 1e-6, 1 - 1e-14):
    for sense in ((PRO,) if spin0 == 0.0 else (PRO, RETRO)):
        a = mp.mpf(spin0) * int(sense)
        x_isco, x_ph = isco(a), photon(a)
        for r in [float(m * x_ph * k) for k in (1.01, 2.0, 40.0)] + [float(m * x_isco)]:
            x = mp.mpf(r) / m
            ut, energy, ang, omega = circular(x, a)
            # 1 - 3/x + 2a/x^(3/2) vanishes at the photon orbit, and nearly so on a near-extremal ISCO
            cond = max(1.0, float(1 / (1 - 3 / x + 2 * a / x ** mp.mpf('1.5'))))
            emit('kerr::circular_orbit_dt_dtau', [KERR_MASS, spin0, sense, r], ut, ULPS * cond)
            emit('kerr::circular_orbit_specific_energy', [KERR_MASS, spin0, sense, r], energy, ULPS * 2 * cond)
            emit('kerr::circular_orbit_specific_angular_momentum', [KERR_MASS, spin0, sense, r], ang * m * c,
                 ULPS * 2 * cond)
            emit('kerr::circular_orbit_angular_velocity', [KERR_MASS, spin0, sense, r], omega * c / m, ULPS)


# ======================= weak field =======================


def shapiro_leg_integral(M, r, b):
    """Integral of 2GM/(c^3 sqrt(s^2 + b^2)) ds from closest approach to the point at distance r."""
    length = mp.sqrt(mp.mpf(r) ** 2 - mp.mpf(b) ** 2)
    if length == 0:
        return mp.mpf(0)
    bb = mp.mpf(b)
    points = [mp.mpf(0)] + [p for p in (bb * 10 ** k for k in range(0, 40, 2)) if p < length] + [length]
    return 2 * g * mp.mpf(M) / c ** 3 * mp.quad(lambda s: 1 / mp.sqrt(s * s + bb * bb), points)


LEGS = [(M_SUN, AU, R_SUN), (M_SUN, 1.523679 * AU, R_SUN), (M_SUN, AU, 5.0 * R_SUN), (M_SUN, 30.07 * AU, 0.3 * AU),
        (M_SUN, AU, AU), (M_SUN, 1.0000001 * AU, AU), (5.9722e24, 3.844e8, 6.4e6), (4.297e6 * M_SUN, 2.5e20, 1e14)]
for M, r, b in LEGS:
    emit('weak_field::shapiro_delay_leg', [M, r, b], shapiro_leg_integral(M, r, b), ULPS)
for M, re, rr, b in [(M_SUN, AU, 1.523679 * AU, R_SUN), (M_SUN, AU, 9.5826 * AU, 1.6 * R_SUN),
                     (M_SUN, AU, 30.07 * AU, 0.1 * AU), (5.9722e24, 2.6561750e7, 4.2164e7, 6.4e6)]:
    emit('weak_field::shapiro_delay', [M, re, rr, b], shapiro_leg_integral(M, re, b) + shapiro_leg_integral(M, rr, b),
         ULPS)


def shapiro_between_integral(M, re, rr, sep):
    """Integral of 2GM/(c^3 |x|) along the straight segment between two points at distances re and rr."""
    re, rr, sep = mp.mpf(re), mp.mpf(rr), mp.mpf(sep)
    # put the emitter at distance re on the x axis; the receiver follows from the cosine rule
    cos_e = (re * re + sep * sep - rr * rr) / (2 * re * sep)     # angle at the emitter between the mass and the path
    foot = re * cos_e                                            # distance along the path to the closest approach
    b2 = re * re - foot * foot
    f = lambda s: 1 / mp.sqrt((s - foot) ** 2 + b2)
    points = sorted({mp.mpf(0), sep} | ({foot} if 0 < foot < sep else set()))
    return 2 * g * mp.mpf(M) / c ** 3 * mp.quad(f, points)


for M, re, rr, sep in [(M_SUN, AU, 1.523679 * AU, 2.2 * AU), (M_SUN, AU, 1.523679 * AU, 0.6 * AU),
                       (M_SUN, AU, AU, 1.0e3), (M_SUN, AU, 5.2038 * AU, 4.2038 * AU), (M_SUN, 2.0 * AU, 3.0 * AU, 4.5 * AU),
                       (5.9722e24, 6.371e6, 2.6561750e7, 2.3e7)]:
    emit('weak_field::shapiro_delay_between', [M, re, rr, sep], shapiro_between_integral(M, re, rr, sep), ULPS)

for M, b in [(M_SUN, R_SUN), (M_SUN, 2.5 * R_SUN), (5.9722e24, 6.371e6), (1e12 * M_SUN, 3e20), (1.898e27, 7.1e7)]:
    bb = mp.mpf(b)
    # integral of the transverse acceleration 2GM b/(c^2 (s^2 + b^2)^(3/2)) over the whole straight path
    bend = 2 * g * mp.mpf(M) / c ** 2 * mp.quad(lambda s: bb / (s * s + bb * bb) ** mp.mpf('1.5'),
                                                [-mp.inf, -bb, 0, bb, mp.inf])
    emit('weak_field::deflection_angle', [M, b], bend, ULPS)

KPC = 3.085_677_581_491_367_3e19
for M, dl, ds, dls in [(M_SUN, 4.0 * KPC, 8.0 * KPC, 4.0 * KPC), (0.3 * M_SUN, 6.5 * KPC, 8.0 * KPC, 1.5 * KPC),
                       (1e12 * M_SUN, 1e6 * KPC, 2.2e6 * KPC, 1.5e6 * KPC), (M_SUN, 550.0 * AU, 1e3 * KPC, 1e3 * KPC)]:
    dl_, ds_, dls_ = mp.mpf(dl), mp.mpf(ds), mp.mpf(dls)
    alpha = lambda theta: 4 * g * mp.mpf(M) / (theta * dl_ * c ** 2)
    # a source exactly behind the lens: theta d_s = alpha(theta d_l) d_ls
    theta = mp.findroot(lambda t: t * ds_ - alpha(t) * dls_, mp.sqrt(4 * g * mp.mpf(M) * dls_ / (c ** 2 * dl_ * ds_)))
    emit('weak_field::einstein_angle', [M, dl, ds, dls], theta, ULPS)

for M, r in [(M_SUN, R_SUN), (5.9722e24, 6.371e6), (5.9722e24, 2.6561750e7), (M_SUN, AU), (1.898e27, 7.1492e7)]:
    u = g * mp.mpf(M) / (mp.mpf(r) * c ** 2)
    emit('weak_field::static_rate_offset', [M, r], -u, ULPS)
    emit('weak_field::circular_orbit_rate_offset', [M, r], -3 * u / 2, ULPS)
for M, r, R, v in [(5.9722e24, 2.6561750e7, 6.378137e6, 465.1), (5.9722e24, 2.6561750e7, 6.371e6, 0.0),
                   (5.9722e24, 6.778e6, 6.378137e6, 465.1), (5.9722e24, 4.2164e7, 6.378137e6, 465.1),
                   (6.4171e23, 2.0428e7, 3.3962e6, 240.7)]:
    um = g * mp.mpf(M) / c ** 2
    exact = um / mp.mpf(R) + mp.mpf(v) ** 2 / (2 * c ** 2) - 3 * um / (2 * mp.mpf(r))
    # conditioning: the three terms partly cancel
    cond = float((um / mp.mpf(R) + mp.mpf(v) ** 2 / (2 * c ** 2) + 3 * um / (2 * mp.mpf(r))) / abs(exact))
    emit('weak_field::orbit_vs_surface_rate_offset', [M, r, R, v], exact, ULPS * max(1.0, cond))

# ======================= Hawking =======================
PHOTON, PAGE, MACGIBBON, COEFFICIENT = 0.0, 1.0, 2.0, 3.0
MACGIBBON_RATE = mp.mpf(5.34e16)


def alpha_of(model, param):
    if model == PHOTON:
        return 1 / (15360 * pi)
    if model == PAGE:
        return mp.mpf(2.011e-4) * (mp.mpf(0.167) + mp.mpf(0.019))
    if model == MACGIBBON:
        return MACGIBBON_RATE * mp.mpf(param) * g ** 2 / (hbar * c ** 4)
    return mp.mpf(param)


HOLES = [1e5, 1e9, 1e12, 7.342e22, M_SUN, 4.297e6 * M_SUN, 1e11 * M_SUN]
MODELS = [(PHOTON, 0.0), (PAGE, 0.0), (MACGIBBON, 1.0), (MACGIBBON, 12.0), (COEFFICIENT, 2.011e-4)]
for M in HOLES:
    mm = mp.mpf(M)
    temp = hbar * c ** 3 / (8 * pi * g * mm * kb)
    area = 4 * pi * (2 * g * mm / c ** 2) ** 2
    emit('hawking::temperature', [M], temp, ULPS)
    emit('hawking::entropy', [M], kb * c ** 3 * area / (4 * g * hbar), ULPS)
    for model, param in MODELS:
        alpha = alpha_of(model, param)
        if model == PHOTON:
            sigma = pi ** 2 * kb ** 4 / (60 * hbar ** 3 * c ** 2)        # Stefan-Boltzmann constant, exact form
            p = sigma * area * temp ** 4
        else:
            p = alpha * hbar * c ** 6 / (g * mm) ** 2
        emit('hawking::power', [M, model, param], p, ULPS)
        emit('hawking::mass_loss_rate', [M, model, param], p / c ** 2, ULPS)
        rate = lambda x: alpha * hbar * c ** 4 / (g * x) ** 2           # -dM/dt at mass x
        life = mp.quad(lambda x: 1 / rate(x), [0, mm])
        emit('hawking::lifetime', [M, model, param], life, ULPS)
for t in [1.0, JULIAN_YEAR, 1e10 * JULIAN_YEAR, 4.35e17, 1e100 * JULIAN_YEAR]:
    for model, param in MODELS:
        alpha = alpha_of(model, param)
        mass = mp.cbrt(3 * alpha * hbar * c ** 4 * mp.mpf(t) / g ** 2)
        # the inverse of the lifetime integral above: check it by integrating back
        back = mp.quad(lambda x: (g * x) ** 2 / (alpha * hbar * c ** 4), [0, mass])
        assert abs(back / mp.mpf(t) - 1) < mp.mpf('1e-100')
        emit('hawking::mass_for_lifetime', [t, model, param], mass, ULPS)
for model, param in MODELS:
    emit('hawking::power_coefficient', [model, param], alpha_of(model, param), ULPS)

# ======================= Newtonian =======================
BODIES = [(3.3011e23, 2.4397e6), (5.9722e24, 6.3710e6), (1.89813e27, 6.9911e7), (M_SUN, R_SUN), (7.342e22, 1.7374e6)]
for M, R in BODIES:
    mm, rr = mp.mpf(M), mp.mpf(R)
    rho = mm / (mp.mpf(4) / 3 * pi * rr ** 3)
    # assemble the sphere shell by shell: dU = G m(r) dm / r
    binding = mp.quad(lambda r: g * (mp.mpf(4) / 3 * pi * r ** 3 * rho) * (4 * pi * r * r * rho) / r, [0, rr])
    emit('newtonian::binding_energy_uniform_sphere', [M, R], binding, ULPS)
    v = mp.sqrt(g * mm / rr)
    emit('newtonian::circular_speed', [M, R], v, ULPS)
    emit('newtonian::escape_speed', [M, R], mp.sqrt(2) * v, ULPS)
    emit('newtonian::circular_period', [M, R], 2 * pi * rr / v, ULPS)
for a_au in [0.387098, 1.0, 1.523679, 5.2038, 30.07]:
    r = a_au * AU
    v = mp.sqrt(g * mp.mpf(M_SUN) / mp.mpf(r))
    emit('newtonian::circular_speed', [M_SUN, r], v, ULPS)
    emit('newtonian::circular_period', [M_SUN, r], 2 * pi * mp.mpf(r) / v, ULPS)
for pa, pb in [(365.256363, 686.980), (686.980, 365.256363), (87.969, 224.701), (4332.589, 10759.22),
               (86400.0, 86164.0905), (1.0, 1.0 + 1e-6)]:
    exact = 1 / abs(1 / mp.mpf(pa) - 1 / mp.mpf(pb))
    emit('newtonian::synodic_period', [pa, pb], exact, ULPS)
for M, r1, r2 in [(M_SUN, AU, 1.523679 * AU), (M_SUN, AU, 5.2038 * AU), (M_SUN, AU, 30.07 * AU),
                  (M_SUN, 1.523679 * AU, AU), (M_SUN, AU, 0.387098 * AU), (5.9722e24, 6.678e6, 4.2164e7),
                  (5.9722e24, 6.778e6, 6.779e6), (M_SUN, AU, AU)]:
    mu, a1, a2 = g * mp.mpf(M), mp.mpf(r1), mp.mpf(r2)
    sma = (a1 + a2) / 2
    emit('newtonian::hohmann_transfer_time', [M, r1, r2], pi * mp.sqrt(sma ** 3 / mu), ULPS)
    vis_viva = lambda r: mp.sqrt(mu * (2 / r - 1 / sma))
    emit('newtonian::hohmann_departure', [M, r1, r2], vis_viva(a1) - mp.sqrt(mu / a1), ULPS)
    emit('newtonian::hohmann_arrival', [M, r1, r2], mp.sqrt(mu / a2) - vis_viva(a2), ULPS)

# ---- write ----
header = ['# Generated by tools/crosscheck.py with mpmath at 120 digits. Do not edit by hand.',
          'function,x1,x2,x3,x4,expected,rel_tol']
constants = [f'constant:{name},{value!r},,,,,' for name, value in CONSTANTS]
text = '\n'.join(header + constants + rows) + '\n'
if '--check' in sys.argv:
    if OUT.read_text(encoding='utf-8') != text:
        print(f'{OUT} is out of date: rerun tools/crosscheck.py')
        sys.exit(1)
    print(f'{OUT} is up to date ({len(rows)} rows)')
else:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding='utf-8')
    print(f'wrote {len(rows)} rows to {OUT}')
