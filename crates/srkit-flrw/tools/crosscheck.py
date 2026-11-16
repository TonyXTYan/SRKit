"""Independent reference values for srkit-flrw.

Writes `tests/golden.csv`, which `tests/golden.rs` reads. Every reference value is computed with mpmath at 50
digits, by quadrature of the Friedmann integrals wherever the Rust crate uses a hypergeometric closed form, so
the two share no code path. Before writing anything the script checks that

  * the hypergeometric closed forms the crate documents agree with that quadrature to 1e-30,
  * mpmath and astropy (`FlatLambdaCDM` with radiation and neutrinos switched off) agree on ages, lookback
    times, Hubble rates and comoving distances.

Run from anywhere:

    python tools/crosscheck.py            # regenerate tests/golden.csv
    python tools/crosscheck.py --check    # exit 1 if the checked-in file differs from a fresh run

Needs mpmath; astropy is optional (its rows and its comparison are skipped without it, and `--check` then
compares the mpmath rows only).

Inputs are written with `repr(float)`, the shortest decimal that round-trips, and converted back to mpmath
from the binary64 value, so Rust and mpmath start from bit-identical numbers. All cosmologies are in Julian
years and light-years.
"""
import math
import sys
from functools import lru_cache
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

GOLDEN = Path(__file__).resolve().parent.parent / "tests" / "golden.csv"

# SI constants, the same values as srkit_core::units
C = 299_792_458.0
JULIAN_YEAR = 31_557_600.0
LIGHT_YEAR = C * JULIAN_YEAR
PARSEC = 3.085_677_581_491_367_3e16
G_NEWTON = 6.674_30e-11
HBAR = 1.054_571_817_646_156_5e-34
BOLTZMANN = 1.380_649e-23
G0 = 9.806_65
KM_S_MPC = 1.0e3 / (1.0e6 * PARSEC)

THIRD = mp.mpf(1) / 3


class Universe:
    """Flat matter + Lambda cosmology in mpmath. Times in years, lengths in light-years."""

    def __init__(self, h0, omega_m):
        self.h0_float, self.om_float = float(h0), float(omega_m)
        self.h0 = mp.mpf(self.h0_float)
        self.om = mp.mpf(self.om_float)
        # the crate forms Omega_Lambda as the binary64 difference 1 - Omega_m
        self.ol = mp.mpf(1.0 - self.om_float)
        self.hl = self.h0 * mp.sqrt(self.ol)
        self.aeq = (self.om / self.ol) ** THIRD

    # ---- elementary closed forms ----
    def a(self, t):
        return self.aeq * mp.sinh(1.5 * self.hl * t) ** (2 * THIRD)

    def ln_a(self, t):
        return mp.log(self.a(t))

    def hubble(self, t):
        return self.hl / mp.tanh(1.5 * self.hl * t)

    def hubble_a(self, a):
        return self.h0 * mp.sqrt(self.om / a**3 + self.ol)

    def t_closed(self, a):
        return mp.asinh((a / self.aeq) ** 1.5) / (1.5 * self.hl)

    # ---- quadrature of the Friedmann integrals (a = s^2 removes the square-root singularity at 0) ----
    def _root(self, s):
        return self.h0 * mp.sqrt(self.om + self.ol * s**6)

    @staticmethod
    def _breaks(upper):
        pts = [mp.mpf(0)]
        p = mp.mpf(1)
        while p < upper:
            pts.append(p)
            p *= 4
        pts.append(upper)
        return pts

    @lru_cache(maxsize=None)
    def t_quad(self, a):
        """t(a) = int_0^a da / (a H)."""
        return mp.quad(lambda s: 2 * s**2 / self._root(s), self._breaks(mp.sqrt(a)))

    @lru_cache(maxsize=None)
    def eta_quad(self, a):
        """eta(a) = int_0^a da / (a^2 H)."""
        return mp.quad(lambda s: 2 / self._root(s), self._breaks(mp.sqrt(a)))

    @lru_cache(maxsize=None)
    def rem_quad(self, a):
        """eta_inf - eta(a) = int_a^inf da / (a^2 H), with u = 1/a.

        mpmath's quad stops on an absolute error estimate, so a late epoch (tiny 1/a) is integrated in the
        rescaled variable u = s / a to keep the relative accuracy."""
        if a > 1:
            return mp.quad(lambda s: 1 / (self.h0 * mp.sqrt(self.om * (s / a) ** 3 + self.ol)), [0, 1]) / a
        return mp.quad(lambda u: 1 / (self.h0 * mp.sqrt(self.om * u**3 + self.ol)), self._breaks(1 / a))

    @lru_cache(maxsize=None)
    def eta_inf_quad(self):
        return mp.quad(lambda s: 2 / self._root(s), [0, 1, 4, 16, 64, mp.inf])

    @lru_cache(maxsize=None)
    def drift_quad(self, a):
        """int_0^t a dt = int_0^a da / H."""
        return mp.quad(lambda s: 2 * s**4 / self._root(s), self._breaks(mp.sqrt(a)))

    def drift_between(self, t0, t1):
        """int_t0^t1 a dt by direct quadrature in t (no cancellation for short intervals)."""
        n = 8
        pts = [t0 + (t1 - t0) * mp.mpf(k) / n for k in range(n + 1)]
        return mp.quad(self.a, pts)

    # ---- the hypergeometric closed forms the crate documents (checked against the quadrature) ----
    def eta_closed(self, a):
        x = self.ol / self.om * a**3
        return 2 / (a * self.hubble_a(a)) * mp.hyp2f1(1, 0.5, mp.mpf(7) / 6, x / (1 + x))

    def rem_closed(self, a):
        x = self.ol / self.om * a**3
        return 1 / (a * self.hubble_a(a)) * mp.hyp2f1(1, 0.5, mp.mpf(4) / 3, 1 / (1 + x))

    def eta_inf_closed(self):
        k = mp.gamma(mp.mpf(1) / 6) * mp.gamma(THIRD) / (3 * mp.sqrt(mp.pi))
        return k / (self.h0 * self.om**THIRD * self.ol ** (mp.mpf(1) / 6))

    def drift_closed_early(self, a):
        x = self.ol / self.om * a**3
        return mp.mpf(2) / 5 * a / self.hubble_a(a) * mp.hyp2f1(1, 0.5, mp.mpf(11) / 6, x / (1 + x))

    def drift_closed_late(self, a):
        x = self.ol / self.om * a**3
        k = -(mp.mpf(2) / 5) * mp.gamma(mp.mpf(11) / 6) * mp.gamma(-THIRD) / mp.sqrt(mp.pi)
        return a / self.hubble_a(a) * mp.hyp2f1(1, 0.5, mp.mpf(2) / 3, 1 / (1 + x)) - k * self.aeq / self.hl

    # ---- derived ----
    def age(self):
        return self.t_quad(mp.mpf(1))

    def ln_a_at_rem(self, target):
        """ln a at which the conformal time left equals `target`: Newton on the quadrature, from the de
        Sitter estimate a = 1 / (H_Lambda * target). d(rem)/d(ln a) = -1 / (a H)."""
        x = mp.log(1 / (self.hl * target))
        for _ in range(60):
            a = mp.exp(x)
            f = self.rem_quad(a) - target
            if abs(f) < abs(target) * mp.mpf("1e-45"):
                return x
            x += f * a * self.hubble_a(a)
        raise ArithmeticError("no convergence")


def self_check(universes):
    """The documented closed forms against quadrature, and the identities the crate relies on."""
    worst = mp.mpf(0)

    def close(x, y):
        nonlocal worst
        e = abs(x / y - 1)
        worst = max(worst, e)
        assert e < mp.mpf("1e-30"), (x, y, e)

    for u in universes:
        close(u.eta_inf_quad(), u.eta_inf_closed())
        for a in ("1e-6", "0.01", "0.3", "0.77", "1", "2.5", "40", "1e5"):
            a = mp.mpf(a)
            close(u.t_quad(a), u.t_closed(a))
            close(u.eta_quad(a), u.eta_closed(a))
            close(u.rem_quad(a), u.rem_closed(a))
            close(u.eta_quad(a) + u.rem_quad(a), u.eta_inf_closed())
            close(u.drift_quad(a), u.drift_closed_early(a))
            if a > mp.mpf("0.01"):  # the late form cancels at tiny a; it is only used for x >= 1
                close(u.drift_quad(a), u.drift_closed_late(a))
        # the particle and event horizons are equal at a = a_eq / 2
        close(u.eta_quad(u.aeq / 2), u.eta_inf_closed() / 2)
    # the de Sitter thrust proper time: both closed forms against quadrature
    for h, big_u, u0, t in ((0.37, 3, 0, 2), (0.37, 1e10, 0, 5), (0.37, 0.5, 7, 3), (0.37, 0, 1e6, 0.001), (1.3, 2, 5, 0.7), (0.9, -2, 3, 4)):
        h, big_u, u0, t = (mp.mpf(v) for v in (h, big_u, u0, t))
        close(de_sitter_tau_quad(h, big_u, u0, t), de_sitter_tau_atanh(h, big_u, u0, t))
        close(de_sitter_tau_quad(h, big_u, u0, t), de_sitter_tau_log(h, big_u, u0, t))
    return worst


def de_sitter_u(h, big_u, u0, t):
    return big_u + (u0 - big_u) * mp.exp(-h * t)


def de_sitter_tau_quad(h, big_u, u0, t):
    pts = [mp.mpf(0)]
    p = min(t, 1 / (h * max(abs(big_u), abs(u0), 1))) / 1000
    while p < t:
        pts.append(p)
        p *= 4
    pts.append(t)
    return mp.quad(lambda s: 1 / mp.sqrt(1 + de_sitter_u(h, big_u, u0, s) ** 2), pts)


def de_sitter_tau_atanh(h, big_u, u0, t):
    gu = mp.sqrt(1 + big_u**2)
    s = 1 - mp.exp(-h * t)
    u1 = de_sitter_u(h, big_u, u0, t)
    rho = gu * s / (mp.sqrt(1 + u1**2) + mp.sqrt(1 + u0**2) * (1 - s))
    return 2 * mp.atanh(rho) / (h * gu)


def de_sitter_tau_log(h, big_u, u0, t):
    gu = mp.sqrt(1 + big_u**2)
    u1 = de_sitter_u(h, big_u, u0, t)
    x1 = (mp.asinh(big_u) + mp.asinh(u1)) / 2
    x0 = (mp.asinh(big_u) + mp.asinh(u0)) / 2
    return (h * t + 2 * (mp.log(mp.cosh(x1)) - mp.log(mp.cosh(x0)))) / (h * gu)


def de_sitter_tau_asymptotic(h, big_u, t):
    """From rest, for h t >> 1: tau = t / gamma_U + C, with C the integral of the transient (the oracle's
    SR-21 method)."""
    gu = mp.sqrt(1 + big_u**2)
    f = lambda s: 1 / mp.sqrt(1 + (big_u * (1 - mp.exp(-h * s))) ** 2) - 1 / gu
    scale = 1 / (h * big_u)
    pts = [mp.mpf(0)] + [scale * mp.mpf(10) ** k for k in range(-2, 60) if scale * mp.mpf(10) ** k < 200 / h] + [200 / h]
    return t / gu + mp.quad(f, pts)


class Rows:
    def __init__(self):
        self.rows = []

    def add(self, source, quantity, u, args, expected, tolerance=1e-12):
        xs = [repr(float(x)) for x in args] + [""] * (4 - len(args))
        self.rows.append(
            (source, quantity, repr(u.h0_float), repr(u.om_float), *xs, mp.nstr(mp.mpf(expected), 25, min_fixed=0, max_fixed=0), "%g" % tolerance)
        )

    def text(self):
        head = [
            "# Reference values for tests/golden.rs. Generated by tools/crosscheck.py; do not edit by hand.",
            "# source: mpmath (50-digit quadrature) or astropy (FlatLambdaCDM, Tcmb0 = 0).",
            "# Cosmologies are in Julian years and light-years; h0 is per year. Temperatures in K, masses in kg,",
            "# accelerations in rapidity per year unless the quantity ends in _si.",
            "source,quantity,h0,omega_m,x1,x2,x3,x4,expected,tolerance",
        ]
        return "\n".join(head + [",".join(r) for r in self.rows]) + "\n"


def fl(x):
    """Round to binary64 and return the exact mpmath value of that double."""
    return mp.mpf(float(x))


def tol(cond):
    """Tolerance for a row whose result amplifies input rounding by `cond`.

    Rows with ordinary conditioning are held to 1e-12 exactly; tests/golden.rs reports those separately.
    Beyond that the tolerance follows the conditioning the crate documents, with a margin of about 10 over
    the rounding of one binary64 operation."""
    cond = abs(float(cond))
    return 1e-12 if cond <= 250 else 4e-15 * cond


def build(with_astropy):
    planck = Universe(67.4 * KM_S_MPC * JULIAN_YEAR, 0.315)
    oracle = Universe(1.0 / 17.53e9 / math.sqrt(1.0 - 0.315), 0.315)
    script = Universe(67.7 / 3.0857e19 * JULIAN_YEAR, 0.31)
    light = Universe(1.0e-10, 0.05)
    heavy = Universe(3.0e-11, 0.95)
    universes = [planck, oracle, script, light, heavy]
    worst = self_check(universes)
    print("closed forms vs 50-digit quadrature: worst relative difference %s" % mp.nstr(worst, 3))

    rows = Rows()
    add = rows.add
    for u in universes:
        age = u.age()
        hl = u.hl
        add("mpmath", "age", u, [], age)
        add("mpmath", "conformal_time_limit", u, [], u.eta_inf_quad())
        add("mpmath", "h_lambda", u, [], hl)
        add("mpmath", "scale_factor_at_equality", u, [], u.aeq)
        add("mpmath", "de_sitter_temperature", u, [], mp.mpf(HBAR) * hl / mp.mpf(JULIAN_YEAR) / (2 * mp.pi * mp.mpf(BOLTZMANN)))
        tds = mp.mpf(HBAR) * hl / mp.mpf(JULIAN_YEAR) / (2 * mp.pi * mp.mpf(BOLTZMANN))
        tcmb = mp.mpf(2.7255)
        add("mpmath", "cmb_crossover_scale_factor", u, [2.7255], tcmb / tds)
        add("mpmath", "cmb_crossover_time", u, [2.7255], u.t_closed(tcmb / tds))
        add("mpmath", "last_light_gain", u, [10.0], mp.log(10) / (3 * hl))
        ln_max = mp.log(mp.mpf(sys.float_info.max))
        add("mpmath", "scale_factor_overflow_time", u, [], u.t_closed(mp.exp(ln_max)))
        r = fl(5.0e4 * LIGHT_YEAR)
        add("mpmath", "lambda_tidal_acceleration_si", u, [r], (hl / mp.mpf(JULIAN_YEAR)) ** 2 * r)
        add("mpmath", "lambda_tidal_acceleration", u, [5.0e4], hl**2 * 5.0e4)
        g = fl(G0 / C * JULIAN_YEAR)
        add("mpmath", "drag_limited_gamma", u, [g], mp.sqrt(1 + (g / hl) ** 2))

        # epochs as multiples of the age, from deep in the matter era to deep in the Lambda era
        for m in (1e-7, 1e-4, 0.01, 0.2, 0.5, 0.74, 1.0, 1.5, 3.0, 8.0, 30.0):
            t = fl(age * m)
            a = u.a(t)
            add("mpmath", "scale_factor", u, [t], a)
            add("mpmath", "ln_scale_factor", u, [t], mp.log(a))
            add("mpmath", "hubble", u, [t], u.hubble(t))
            add("mpmath", "conformal_time", u, [t], u.eta_quad(a))
            add("mpmath", "conformal_time_remaining", u, [t], u.rem_quad(a))
            add("mpmath", "particle_horizon_proper", u, [t], a * u.eta_quad(a))
            add("mpmath", "event_horizon_proper", u, [t], a * u.rem_quad(a))
            add("mpmath", "reachable_share_of_observable", u, [t], (u.rem_quad(a) / u.eta_quad(a)) ** 3)
            add("mpmath", "reach_volume_ratio", u, [t, fl(age)], (u.rem_quad(a) / u.rem_quad(u.a(fl(age)))) ** 3)
            f = u.aeq * hl * u.rem_quad(a)
            add("mpmath", "reach_mass", u, [t], mp.mpf(C) ** 3 * f**3 / (2 * mp.mpf(G_NEWTON) * hl / mp.mpf(JULIAN_YEAR)))
            add("mpmath", "scale_factor_integral", u, [t], u.drift_quad(a))
            add("mpmath", "cmb_temperature", u, [t, 2.7255], tcmb / a)

            # inverses: the input is the forward value rounded to binary64; the first-order correction below
            # is exact to the square of that rounding (1e-32)
            af = fl(a)
            add("mpmath", "time_at_scale_factor", u, [af], u.t_quad(af))
            add("mpmath", "hubble_at_scale_factor", u, [af], u.hubble_a(af))
            add("mpmath", "conformal_time_at_scale_factor", u, [af], u.eta_quad(af))
            lf = fl(mp.log(a))
            add("mpmath", "time_at_ln_scale_factor", u, [lf], u.t_quad(mp.exp(lf)))
            e = u.eta_quad(a)
            rm = u.rem_quad(a)
            ef = fl(e)
            # near eta_inf the inverse amplifies the rounding of eta by eta / (eta_inf - eta); past 1e6 the
            # question is no longer meaningful in binary64 and the crate says to use the remaining time
            cond = e / rm
            if cond < 1e6:
                add("mpmath", "time_at_conformal_time", u, [ef], t + (ef - e) * a, tol(cond))
            rf = fl(rm)
            # early on it is the remaining time that sits within rounding of eta_inf: d ln t / d ln rem
            cond = rm * a / t
            add("mpmath", "time_at_conformal_time_remaining", u, [rf], t - (rf - rm) * a, tol(cond))
            add("mpmath", "last_contact_time", u, [rf], t - (rf - rm) * a, tol(cond))
            add("mpmath", "last_reply_time", u, [fl(rf / 2)], t - (2 * fl(rf / 2) - rm) * a, tol(cond))

        # far future, where a itself overflows binary64: everything stays finite through ln a
        for ht in (100.0, 600.0, 5000.0):
            t = fl(ht / hl)
            theta = 1.5 * hl * t
            ln_a = mp.log(u.aeq) + 2 * THIRD * mp.log(mp.sinh(theta))
            add("mpmath", "ln_scale_factor", u, [t], ln_a)
            add("mpmath", "hubble", u, [t], u.hubble(t))
            add("mpmath", "time_at_ln_scale_factor", u, [fl(ln_a)], u.t_closed(mp.exp(fl(ln_a))))
            a = mp.exp(ln_a)
            rm = u.rem_quad(a)
            add("mpmath", "event_horizon_proper", u, [t], a * rm)
            if ht < 700:  # beyond that the conformal time left underflows binary64
                # these fall as exp(-H t), so the rounding of H t is amplified by H t
                add("mpmath", "conformal_time_remaining", u, [t], rm, tol(ht))
                add("mpmath", "cmb_temperature", u, [t, 2.7255], tcmb * mp.exp(-ln_a), tol(ht))
                rf = fl(rm)
                add("mpmath", "last_contact_time", u, [rf], t - (rf - rm) * a)

        # redshifts
        for z in (0.01, 0.5, 1.0, 3.0, 20.0, 1100.0, -0.5, -0.99):
            a = 1 / (1 + mp.mpf(z))
            tz = u.t_quad(a)
            add("mpmath", "time_at_redshift", u, [z], tz)
            # both are differences from today, so a small |z| cancels
            add("mpmath", "lookback_time", u, [z], age - tz, tol(age / (age - tz)))
            if z > -0.9:
                d = u.eta_quad(mp.mpf(1)) - u.eta_quad(a)
                add("mpmath", "comoving_distance_at_redshift", u, [z], d, tol(u.eta_quad(mp.mpf(1)) / d))

        # pairs of epochs
        for m1, m2 in ((0.001, 0.3), (0.3, 0.6), (0.2, 1.0), (1.0, 1.0001), (1.0, 2.0), (1.0, 12.0), (3.0, 9.0), (9.0, 3.0), (0.6, 40.0)):
            t1, t2 = fl(age * m1), fl(age * m2)
            a1, a2 = u.a(t1), u.a(t2)
            # cancellation for close epochs is documented: eta / |d eta|, or the same with the conformal
            # time left when both epochs are after equality
            if a1 > u.aeq and a2 > u.aeq:
                d = u.rem_quad(a1) - u.rem_quad(a2)
                cond = max(u.rem_quad(a1), u.rem_quad(a2)) / d
            else:
                d = u.eta_quad(a2) - u.eta_quad(a1)
                cond = max(u.eta_quad(a1), u.eta_quad(a2)) / d
            add("mpmath", "comoving_distance", u, [t1, t2], d, tol(cond))
            z = a2 / a1 - 1
            add("mpmath", "redshift_between", u, [t1, t2], z, tol((1 + abs(mp.log(a1 / u.aeq)) + abs(mp.log(a2 / u.aeq))) * (1 + 1 / abs(z))))
            for p in (1e-3, 0.7, 30.0, 1e9):
                add("mpmath", "coasting_proper_velocity", u, [p, t1, t2], p * a1 / a2)
                add("mpmath", "coasting_rapidity", u, [mp.asinh(p), t1, t2], mp.asinh(mp.sinh(fl(mp.asinh(p))) * a1 / a2))
                gam = fl(mp.sqrt(1 + mp.mpf(p) ** 2))
                # gamma - 1 is what carries the information when gamma is close to 1
                g2 = mp.sqrt(1 + (gam**2 - 1) * (a1 / a2) ** 2)
                add("mpmath", "coasting_gamma", u, [gam, t1, t2], g2)
            add("mpmath", "coasting_rapidity", u, [300.0, t1, t2], mp.asinh(mp.sinh(300) * a1 / a2))

        # light cones: emit at t1, comoving distance a given fraction of what is left
        for m1 in (0.01, 0.5, 1.0, 4.0, 20.0):
            t1 = fl(age * m1)
            a1 = u.a(t1)
            left = u.rem_quad(a1)
            for frac in (1e-6, 0.01, 0.5, 0.9, 0.999999):
                chi = fl(left * frac)
                rem_r = left - chi
                # receive where the conformal time left is rem_r
                a_r = mp.exp(u.ln_a_at_rem(rem_r))
                t_r = u.t_quad(a_r)
                # conditioning: the rounding of (left - chi), relative to rem_r, times d ln t / d ln rem
                cond = left / rem_r * max(1, rem_r * a_r / t_r)
                add("mpmath", "light_arrival_time", u, [t1, chi], t_r, tol(cond))
                # and back: light received at t_r (rounded) from chi away left at t1, to first order in the
                # rounding of t_r. Conditioning: rounding of the conformal time left, times d ln t / d ln rem
                trf = fl(t_r)
                t_e = t1 - (u.rem_quad(u.a(trf)) - rem_r) * a1
                add("mpmath", "light_emission_time", u, [trf, chi], t_e, tol(left * a1 / t1))
            # the conformal-chart form, for emission before the Lambda era crowds eta against eta_inf
            if m1 <= 1.0:
                e1 = fl(u.eta_quad(a1))
                for frac in (0.01, 0.5, 0.9):
                    chi = fl(left * frac)
                    er = fl(float(e1) + float(chi))  # the sum the crate forms, in binary64
                    rem_r = u.eta_inf_quad() - er
                    add("mpmath", "light_reception", u, [e1, chi], u.t_quad(mp.exp(u.ln_a_at_rem(rem_r))), tol(u.eta_inf_quad() / rem_r))

        # thrust in full Lambda-CDM: u(t1) = [a0 u0 + alpha * int a dt] / a1
        for m0, m1, u0 in ((0.5, 1.0, 0.0), (1.0, 1.01, 0.0), (1.0, 2.0, 0.0), (1.0, 2.0, 5.0), (1.0, 12.0, 1e11), (0.01, 30.0, 0.0), (2.0, 2.5, -3.0)):
            t0, t1 = fl(age * m0), fl(age * m1)
            a0, a1 = u.a(t0), u.a(t1)
            for alpha in (g, fl(g * 1e-9)):
                val = (a0 * u0 + alpha * u.drift_between(t0, t1)) / a1
                # the thrust term is a difference of two integrals from t = 0: documented loss 1e-16 / (H dt)
                add("mpmath", "thrust_proper_velocity", u, [u0, alpha, t0, t1], val, tol(u.drift_quad(a1) / u.drift_between(t0, t1)))

        # thrust in the de Sitter limit
        for u0, alpha, hdt in ((0.0, g, 1e-12), (0.0, g, 1e-6), (0.0, g, 0.3), (0.0, g, 3.0), (0.0, g, 50.0), (0.0, g, 2000.0),
                               (2.0, fl(hl), 0.5), (40.0, fl(hl * 3), 5.0), (1e8, 0.0, 0.01), (1e8, 0.0, 30.0), (5.0, fl(-hl * 2), 1.5),
                               (0.0, fl(hl * 1e-3), 0.7), (1e-3, fl(hl * 1e-3), 1e-5)):
            dt = fl(hdt / hl)
            big_u = mp.mpf(alpha) / hl
            add("mpmath", "de_sitter_thrust_proper_velocity", u, [u0, alpha, dt], de_sitter_u(hl, big_u, mp.mpf(u0), dt))
            if hdt <= 60:
                tau = de_sitter_tau_quad(hl, big_u, mp.mpf(u0), dt)
            else:
                assert u0 == 0.0
                tau = de_sitter_tau_asymptotic(hl, big_u, dt)
            add("mpmath", "de_sitter_thrust_proper_time", u, [u0, alpha, dt], tau)
        # the oracle's SR-21 case: 1 g from rest for 1e100 years
        add("mpmath", "de_sitter_thrust_proper_time", u, [0.0, g, 1e100], de_sitter_tau_asymptotic(hl, g / hl, mp.mpf(1e100)))

    if with_astropy:
        astropy_rows(rows, planck, script)
    return rows


def astropy_rows(rows, *universes):
    """Ages, lookback times, Hubble rates and comoving distances from astropy, radiation switched off."""
    from astropy import units as au
    from astropy.cosmology import FlatLambdaCDM

    worst = {}

    def note(name, value, exact):
        worst[name] = max(worst.get(name, 0.0), abs(value / float(exact) - 1))

    for u in universes:
        h0_km_s_mpc = u.h0_float / JULIAN_YEAR / KM_S_MPC
        cosmo = FlatLambdaCDM(H0=h0_km_s_mpc, Om0=u.om_float, Tcmb0=0.0)
        mpc_ly = 1.0e6 * PARSEC / LIGHT_YEAR
        age = cosmo.age(0).to_value(au.yr)
        rows.add("astropy", "age", u, [], age, 1e-12)
        note("age", age, u.age())
        for z in (0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 1100.0):
            tz = cosmo.age(z).to_value(au.yr)
            lb = cosmo.lookback_time(z).to_value(au.yr)
            dc = cosmo.comoving_distance(z).to_value(au.Mpc) * mpc_ly
            hz = cosmo.H(z).to_value(au.km / au.s / au.Mpc) * KM_S_MPC * JULIAN_YEAR
            a = 1 / (1 + mp.mpf(z))
            ref = {
                "time_at_redshift": (tz, u.t_quad(a), 1e-12),
                "lookback_time": (lb, u.age() - u.t_quad(a), tol(4 / z)),
                "comoving_distance_at_redshift": (dc, u.eta_quad(mp.mpf(1)) - u.eta_quad(a), tol(4 / z)),
                "hubble_at_scale_factor": (hz, u.hubble_a(mp.mpf(float(a))), 1e-12),
            }
            for name, (value, exact, tolerance) in ref.items():
                arg = float(a) if name == "hubble_at_scale_factor" else z
                rows.add("astropy", name, u, [arg], value, tolerance)
                note(name, value, exact)
            # the scale factor at astropy's age: a(t(z)) = 1 / (1 + z)
            rows.add("astropy", "scale_factor", u, [tz], 1 / (1 + z), 1e-12)
    print("astropy (Tcmb0 = 0) vs mpmath, worst relative difference:")
    for name, value in worst.items():
        print("  %-32s %.2e" % (name, value))


def main():
    check = "--check" in sys.argv[1:]
    try:
        import astropy  # noqa: F401

        with_astropy = True
    except ImportError:
        with_astropy = False
        print("astropy not installed: its rows are skipped")
    text = build(with_astropy).text()
    n = text.count("\n") - 5
    if check:
        old = GOLDEN.read_text(encoding="utf-8")
        if not with_astropy:
            old = "".join(line + "\n" for line in old.splitlines() if not line.startswith("astropy,"))
        if old != text:
            print("tests/golden.csv differs from a fresh run")
            sys.exit(1)
        print("tests/golden.csv matches a fresh run (%d rows)" % n)
    else:
        GOLDEN.write_text(text, encoding="utf-8")
        print("wrote %s (%d rows)" % (GOLDEN, n))


if __name__ == "__main__":
    main()
