"""The physics oracle (`tools/physics-checks/oracle.csv`): every row that is special relativity, optics or exact
timekeeping, recomputed through the Python API and compared at the row's own tolerance.

These are the rows the Rust core's `tests/oracle.rs` reproduces, computed the same way. The cosmology and
gravity rows belong to other crates, and the game-design rows to `oracle_check.py`.
"""

from __future__ import annotations

import math

import pytest

import srkit
from helpers import D_PROXIMA, G, INF, ORIGIN, X_AXIS
from srkit import causal, kinematics, optics, rocket
from srkit.units import AU, C, G0, JULIAN_YEAR, LIGHT_YEAR, Units

J_PER_KT = 4.184e12

#: The rows computed below, in the order of the table.
ROWS = [
    *(f"SR-{n:02d}" for n in range(1, 20)),
    "SR-22",
    *(f"SR-{n}" for n in range(23, 33)),
    "SR-50",
    "SR-51",
    "SR-52",
    "SR-63",
    "SR-69",
    "SR-70",
    "SR-72",
    "SR-73",
    "SR-80",
    "SR-81",
    "SR-82",
    "SR-87",
    "DT-02",
    "DT-11",
    "DT-12",
    "DT-13",
    "DT-14",
]


def causal_failures(samples: int = 100_000) -> tuple[int, int]:
    """The oracle's DT-14 experiment: a receiver moving at up to 0.9 c, placed so that the signal arrives a hair
    after a whole tick. `ceil` of the floating-point closed form is sometimes a tick early; settling it with
    `causal.settle_arrival` must never be. Returns `(failures, times the plain ceiling was early)`."""
    c_m_per_ms = 299_792.458
    mask = (1 << 64) - 1
    state = 0x2545F4914F6CDD1D  # xorshift, fixed seed

    def uniform() -> float:
        nonlocal state
        state ^= (state << 13) & mask
        state ^= state >> 7
        state ^= (state << 17) & mask
        return (state >> 11) / (1 << 53)

    failures = naive_early = 0
    emission = srkit.Time384.from_decimal("1" + "0" * 39)
    for _ in range(samples):
        k = math.floor(1e6 + uniform() * (1e13 - 1e6))
        v = (uniform() * 1.8 - 0.9) * c_m_per_ms  # receiver velocity along the line, m per ms
        target = k + uniform() * 4.0 * (k * 2.0**-52)
        x0 = target * (c_m_per_ms - v)  # receiver's offset at emission, m

        def displacement(tick, x0=x0, v=v):
            dt = tick.local_seconds_since(emission) * 1e3
            return causal.grid_from_metres(x0 + v * dt, 0.0, 0.0)

        estimate = emission + math.ceil(x0 / (c_m_per_ms - v))
        if not causal.arrived(emission, estimate, displacement(estimate)):
            naive_early += 1
        settled = causal.settle_arrival(emission, estimate, displacement, guard=8)
        on_time = causal.arrived(emission, settled, displacement(settled))
        earlier = settled - 1
        minimal = not causal.arrived(emission, earlier, displacement(earlier))
        if not (on_time and minimal):
            failures += 1
    return failures, naive_early


def compute() -> dict[str, float]:
    yr = Units.YEARS
    rows: dict[str, float] = {}

    # the scales of a 10 m/s^2 drive
    rows["SR-01"] = 1.0 / yr.accel_from_si(10.0)
    rows["SR-02"] = yr.length_from_si(C * C / 10.0)

    # 1 g flip-and-burn across the galaxy
    trip = srkit.FlipBurn(1e5, G)
    rows["SR-03"] = trip.tau
    rows["SR-04"] = trip.t
    rows["SR-05"] = trip.gamma_peak
    rows["SR-06"] = trip.mass_ratio(1.0)
    rows["SR-07"] = 2.0 * trip.phi_peak
    # the same ride as a plan of legs must agree with the closed form
    legs = srkit.Trip(1e5, G)
    assert legs.tau == pytest.approx(trip.tau, rel=1e-14)
    assert legs.lag == pytest.approx(trip.lag, rel=1e-12)

    # proper time to reach 0.9 c at 1 g: evaluate the burn where its speed is 0.9 c
    burn = srkit.Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, G)
    t_09 = 0.9 / math.sqrt(1.0 - 0.81) / G  # gamma beta / a
    assert burn.at(t_09)[1][0] == pytest.approx(0.9, abs=1e-15)
    rows["SR-08"] = burn.dtau(t_09)

    # loops flown at 1 g
    for row_id, r in [("SR-09", 1.0), ("SR-10", 100.0), ("SR-11", 1e4), ("SR-12", 5e4)]:
        rows[row_id] = srkit.CircularOrbit.from_radius_acceleration(r, G).gamma
    ring = srkit.CircularOrbit.from_radius_acceleration(100.0, G)
    rows["SR-13"] = ring.period
    rows["SR-14"] = ring.lap_tau
    rows["SR-15"], rows["SR-16"] = ring.doppler_range()
    # lap average of the Doppler factor of a distant source in the plane, over gamma
    n = 3600
    angles = (math.tau * (k + 0.5) / n for k in range(n))
    mean = sum(optics.observe(ring.phi, math.cos(th), math.sin(th), 0.0).doppler for th in angles) / n
    rows["SR-17"] = mean / ring.gamma
    rows["SR-18"] = ring.sagnac()
    rows["SR-19"] = ring.sagnac() * ring.gamma

    # 1 g held for 1e100 years of coordinate time, in flat space
    rows["SR-22"] = burn.dtau(1e100)
    rows["DT-12"] = burn.dtau(1e100) * G

    # warning time: light from a launch 20 light-minutes away against the projectile itself
    minutes = Units(60.0)
    d = 20.0
    for row_id, beta, unit_s in [
        ("SR-23", 0.01, 3600.0),
        ("SR-24", 0.5, 60.0),
        ("SR-25", 0.9999, 1.0),
        ("SR-26", 0.99999, 1e-3),
    ]:
        shot = srkit.Segment.inertial(0.0, INF, (d, 0.0, 0.0), (-beta, 0.0, 0.0))
        target = srkit.Segment.rest(-INF, INF, ORIGIN)
        seen = target.recv(0.0, (d, 0.0, 0.0))
        # when the shot reaches the target: its light, sent from the target's position, "arrives" at once
        hit = shot.send(d / beta, ORIGIN, 0.0, INF)
        assert hit == pytest.approx(d / beta, rel=1e-9)
        # the shot covers the last (1 - beta) of the light's lead: d / beta - d = d (1 - beta) / beta
        warning = d * kinematics.speed_shortfall(kinematics.rapidity_of_speed(beta)) / beta
        assert abs(warning - (hit - seen)) <= 1e-9 * hit
        rows[row_id] = minutes.time_to_si(warning) / unit_s

    # kinetic energy, gamma and head-on blueshift near c
    for beta, ke, gam, blue in [(0.9999, "SR-27", "SR-28", "SR-29"), (0.99999, "SR-30", "SR-31", "SR-32")]:
        phi = kinematics.rapidity_of_speed(beta)
        rows[ke] = rocket.kinetic_per_rest_energy(phi) * C * C / (J_PER_KT * 1e6)
        rows[gam] = kinematics.gamma_of_speed(beta)
        rows[blue] = kinematics.doppler_of_rapidity(phi)
    for row_id, beta in [("SR-50", 0.01), ("SR-51", 0.1)]:
        phi = kinematics.rapidity_of_speed(beta)
        rows[row_id] = rocket.kinetic_per_rest_energy(phi) * C * C / J_PER_KT

    # a slow ship overtaken by a later, faster one: the gap between their arrivals at Proxima
    slow = srkit.Segment.inertial(0.0, INF, ORIGIN, (0.01, 0.0, 0.0))
    fast = srkit.Segment.inertial(0.0, INF, ORIGIN, (0.1, 0.0, 0.0))

    def arrival(seg: srkit.Segment) -> float:
        return D_PROXIMA / seg.at(0.0)[1][0]

    assert slow.at(arrival(slow))[0][0] == pytest.approx(D_PROXIMA, abs=1e-12)
    rows["SR-52"] = arrival(slow) - arrival(fast)
    spun_up = srkit.CircularOrbit.from_radius_acceleration(5e4, G)
    rows["SR-63"] = (round(spun_up.gamma) - 1.0) * 1e6 * C * C

    # beams and apertures
    rows["SR-69"] = optics.rayleigh_resolution(1e6, 1e-6, 100.0 * LIGHT_YEAR) / 1.22 / 1e3
    rows["SR-70"] = 1e6 * 1e6 / 1e-6 / LIGHT_YEAR
    rows["SR-87"] = optics.rayleigh_resolution(AU, 1e-6, 1.0) / 1.22

    # the galaxy's rotation
    beta_rot = Units.SECONDS.speed_from_si(220e3)
    rows["SR-72"] = beta_rot
    rows["SR-73"] = rocket.kinetic_per_rest_energy(kinematics.rapidity_of_speed(beta_rot))

    # a 1,000 t ship to Proxima at 1 g: a beamed sail to the turnover, and a photon rocket for the whole ride
    prox = srkit.FlipBurn(D_PROXIMA, G)
    rest_energy = 1e6 * C * C
    # a sail held at 1 g receives its beam redshifted twice over (D^2, D = e^-phi): E = (m c^2 / 2)(e^phi - 1)
    sail_energy = 0.5 * rest_energy * math.expm1(prox.phi_peak)
    assert rocket.sail_rapidity(sail_energy / rest_energy) == pytest.approx(prox.phi_peak, abs=1e-12)
    rows["SR-80"] = sail_energy
    rows["SR-81"] = 0.5 * 1e6 * G0 * C * math.exp(2.0 * prox.phi_peak)
    rows["SR-82"] = rocket.propellant_per_dry_mass(2.0 * prox.phi_peak, 1.0) * rest_energy

    # playback rate of the launch site's clock at rapidity 5 on the way out
    far = srkit.Trip(1e6, G)
    rows["DT-11"] = far.state(5.0 / G).doppler_home

    # exact clocks: ranges of the stored integer, in Julian years
    ms_per_year = JULIAN_YEAR * 1e3
    rows["DT-02"] = int(srkit.Time384.MAX) / ms_per_year
    assert int(srkit.Time512.MAX) / ms_per_year > 1e140
    # an orbit's phase at 1e20 years by integer modular reduction is exact: the same phase as at the remainder
    period_ms = 31_557_600_000
    late = srkit.Time384.from_decimal("3155760000000" + "0" * 17)
    phase = (late + 123_456_789) % period_ms
    rows["DT-13"] = (phase - 123_456_789) / period_ms

    # causal rounding: adversarial arrivals just above a tick, settled by the exact integer check
    failures, naive_early = causal_failures()
    assert naive_early > 0, "the experiment should catch the plain ceiling out at least once"
    rows["DT-14"] = float(failures)
    return rows


@pytest.fixture(scope="module")
def computed() -> dict[str, float]:
    return compute()


@pytest.mark.parametrize("row_id", ROWS)
def test_oracle_row(row_id, oracle, computed):
    row = oracle[row_id]
    value = computed[row_id]
    assert row.agrees(value), f"{row!r}; computed {value!r}"


def test_the_rows_computed_are_the_rows_listed(computed):
    assert sorted(computed) == sorted(ROWS)
    assert len(ROWS) == 47  # the same 47 rows as the Rust core's oracle test


def test_the_named_rows_to_full_precision(computed):
    """The closed forms behind a few rows, to rounding rather than to the table's two or three figures."""
    a, d = G, 1e5
    assert computed["SR-03"] == pytest.approx(2.0 / a * math.acosh(1.0 + a * d / 2.0), rel=1e-14)
    assert computed["SR-04"] == pytest.approx(math.sqrt(d * d + 4.0 * d / a), rel=1e-15)
    assert computed["SR-05"] == pytest.approx(1.0 + a * d / 2.0, rel=1e-15)
    assert computed["SR-08"] == pytest.approx(math.atanh(0.9) / a, rel=1e-14)
    assert computed["SR-22"] == pytest.approx(math.asinh(a * 1e100) / a, rel=1e-14)
    assert computed["DT-12"] == pytest.approx(math.log(2.0 * a * 1e100), rel=1e-14)
    assert computed["SR-28"] == pytest.approx(1.0 / math.sqrt(1.0 - 0.9999**2), rel=1e-12)
    assert computed["SR-29"] == pytest.approx(math.sqrt(1.9999 / 0.0001), rel=1e-12)
    assert computed["SR-31"] == pytest.approx(1.0 / math.sqrt(1.0 - 0.99999**2), rel=1e-11)
    assert computed["SR-32"] == pytest.approx(math.sqrt(1.99999 / 0.00001), rel=1e-11)
    # a loop held at proper acceleration a: gamma^2 beta^2 = a r, so gamma = sqrt(1 + a r)
    for row_id, r in [("SR-09", 1.0), ("SR-10", 100.0), ("SR-11", 1e4), ("SR-12", 5e4)]:
        assert computed[row_id] == pytest.approx(math.sqrt(1.0 + a * r), rel=1e-14)
