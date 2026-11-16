"""Itineraries, rides in null coordinates, orbits, and the function modules."""

from __future__ import annotations

import math
import random

import pytest

import srkit
from helpers import D_PROXIMA, G, INF, ORIGIN, X_AXIS, point, scaled, unit
from srkit import Segment, Worldline, kepler, kinematics, optics, plan, rocket


# ------------------------------------------------------------------------------------------------- plan


def test_a_line_joins_and_arrives():
    """A straight run reaches its target, joins continuously and spends the rapidity it reports."""
    rng = random.Random(41)
    planned = refused = 0
    for _ in range(300):
        x1 = (rng.uniform(-50, 50), rng.uniform(-50, 50), rng.uniform(-5, 5))
        if math.hypot(*x1) <= 0.5:
            continue
        a, brake_ratio = rng.uniform(0.05, 3.0), rng.uniform(0.3, 3.0)
        cap = INF if rng.random() < 0.5 else rng.uniform(0.2, 3.0)
        phi0 = 0.0 if rng.random() < 0.5 else rng.uniform(0.0, 0.15)
        phi1 = 0.0 if rng.random() < 0.5 else rng.uniform(0.0, 0.15)
        start = plan.Start(t=2000.0, tau=7.0)
        try:
            line = plan.line(start, x1, a, cap, phi0, phi1, a * brake_ratio)
        except ValueError:
            # refused: only when the run is too short for the speeds asked
            assert phi0 > 0.0 or phi1 > 0.0
            refused += 1
            continue
        planned += 1
        length = math.hypot(*x1)
        assert line.pieces and len(line.pieces) == len(line.segments)
        assert [p.segment for p in line.pieces] == list(line.segments)
        Worldline(line.segments, anchor_t=2000.0, anchor_tau=7.0).check_continuity(1e-9 * (1.0 + length), 1e-9)
        last = line.pieces[-1].end()
        assert math.dist(last.x, x1) <= 1e-9 * (1.0 + length)
        assert line.end.x == x1  # pinned on the target
        assert last.phi == pytest.approx(phi1, abs=1e-9)
        assert max(phi0, phi1) - 1e-12 <= line.phi_peak <= cap + 1e-12
        assert line.length == pytest.approx(length, rel=1e-15)
        assert math.hypot(*line.n) == pytest.approx(1.0)
        # proper time adds up, and the phases come in order
        total = sum(p.dtau for p in line.pieces)
        assert line.end.tau - 7.0 == pytest.approx(total, rel=1e-9, abs=1e-9)
        phases = [p.phase for p in line.pieces]
        assert all(a != b for a, b in zip(phases, phases[1:]))
        assert set(phases) <= {"boost", "coast", "brake"}
        assert "coast" not in phases or line.coast_length > 0.0
        assert line.dphi == pytest.approx((line.phi_peak - phi0) + (line.phi_peak - phi1))
    assert planned > 200 and planned + refused == 300


def test_a_line_agrees_with_the_trip_planner():
    """From rest to rest, a line is the flip-and-burn (or the capped leg) of the trip planner."""
    rng = random.Random(42)
    for _ in range(200):
        d = math.exp(rng.uniform(math.log(0.01), math.log(1e5)))
        a = rng.uniform(0.05, 3.0)
        cap = INF if rng.random() < 0.5 else rng.uniform(1.05, 50.0)
        line = plan.line(plan.Start(), (d, 0.0, 0.0), a, math.acosh(cap) if cap < INF else INF)
        trip = srkit.Trip(d, a, cap)
        assert line.end.tau == pytest.approx(trip.tau, rel=1e-9)
        assert line.end.t == pytest.approx(trip.t_end, rel=1e-9)
        assert line.dphi == pytest.approx(2.0 * trip.leg.phc, rel=1e-9, abs=1e-9)
        assert trip.phi_peak == trip.leg.phc and trip.gamma_peak == pytest.approx(math.cosh(trip.phi_peak))
        if cap == INF:
            fb = srkit.FlipBurn(d, a)
            assert trip.tau == pytest.approx(fb.tau, rel=1e-13)
            assert trip.t_end == pytest.approx(fb.t, rel=1e-13)
            assert srkit.flip_burn_tau(d, a) == pytest.approx(fb.tau, rel=1e-13)
            assert srkit.flip_burn_distance(fb.tau, a) == pytest.approx(d, rel=1e-9)


def test_short_hops_keep_their_digits():
    # One astronomical unit at 1 g: the peak rapidity is 4e-3, where acosh(1 + x) would lose half its digits.
    au = srkit.units.AU_LY
    line = plan.line(plan.Start(), (au, 0.0, 0.0), G)
    # Newtonian limit: tau = 2 sqrt(d / a), corrected by 1 - a d / 24 to leading order
    expected = 2.0 * math.sqrt(au / G) * (1.0 - G * au / 24.0)
    assert line.end.tau == pytest.approx(expected, rel=1e-10)
    assert math.dist(line.pieces[1].end().x, (au, 0.0, 0.0)) < 1e-15 * 1e-3


def test_pieces_and_degenerate_plans():
    start = plan.Start(t=5.0, tau=1.0, x=(1.0, 2.0, 3.0))
    assert (start.t, start.tau, start.x) == (5.0, 1.0, (1.0, 2.0, 3.0))
    assert repr(start) == "Start(t=5.0, tau=1.0, x=(1.0, 2.0, 3.0))"
    assert repr(plan.Start()) == "Start(t=0.0, tau=0.0, x=(0.0, 0.0, 0.0))"

    nowhere = plan.line(start, start.x, 1.0)
    assert nowhere.pieces == () and nowhere.segments == () and nowhere.dphi == 0.0
    assert nowhere.end.x == start.x and nowhere.end.t == 5.0

    jump = plan.null_piece(start, (1.0, 2.0, 7.0))
    assert jump.dtau == 0.0 and jump.phase == "jump" and jump.segment.t1 == 9.0
    assert jump.end().x == (1.0, 2.0, 7.0) and jump.end().tau == 1.0

    wait = plan.rest(start, 2.5)
    assert wait.phase == "rest" and wait.tau0 == 1.0 and wait.end().tau == 3.5 and wait.end().t == 7.5
    assert wait.segment == Segment.rest(5.0, 7.5, (1.0, 2.0, 3.0))
    assert plan.rest(start, INF).segment.t1 == INF

    coast = plan.lin_piece(start, (0.6, 0.0, 0.0), 4.0)
    assert coast.phase == "coast" and coast.end().t == pytest.approx(5.0 + 4.0 / 0.8)
    assert coast.end().x == pytest.approx((1.0 + 0.6 * 5.0, 2.0, 3.0))

    boost = plan.hyp_piece(start, X_AXIS, 0.5, 0.0, 1.0)
    assert boost.phase == "boost" and boost.dtau == 2.0 and boost.end().phi == pytest.approx(1.0)
    slow = plan.hyp_piece(boost.end(), X_AXIS, -0.5, 1.0, 0.0, "brake")
    assert slow.phase == "brake" and slow.end().v == pytest.approx((0.0, 0.0, 0.0), abs=1e-15)
    assert "phase='brake'" in repr(slow) and "Segment.hyperbolic(" in repr(slow)
    # an End carries on as the next Start, directly or through .start()
    assert plan.rest(boost.end(), 1.0).segment == plan.rest(boost.end().start(), 1.0).segment
    assert repr(boost.end()).startswith("End(t=")

    stop = plan.brake(start, X_AXIS, 0.0, 1.0)
    assert stop.pieces == () and stop.segments == () and stop.dphi == 0.0 and stop.end.x == start.x
    stop = plan.brake(start, X_AXIS, 1.2, 0.6)
    assert len(stop.pieces) == 1 and stop.dphi == 1.2 and stop.pieces[0].dtau == pytest.approx(2.0)
    assert stop.end.phi == pytest.approx(0.0, abs=1e-15) and stop.segments[0] == stop.pieces[0].segment

    for bad in [
        lambda: plan.line(start, (9.0, 2.0, 3.0), 0.0),
        lambda: plan.line(start, (9.0, 2.0, 3.0), 1.0, phi_cap=0.0),
        lambda: plan.line(start, (9.0, 2.0, 3.0), 1.0, phi0=-0.1),
        lambda: plan.line(start, (9.0, 2.0, 3.0), 1.0, a_brake=-1.0),
        lambda: plan.line(start, (1.0, 2.0, 3.0 + 1e-6), 1.0, phi0=0.5),  # too short to shed that speed
        lambda: plan.hyp_piece(start, X_AXIS, 0.5, 1.0, 0.0),  # the wrong sign of thrust
        lambda: plan.hyp_piece(start, X_AXIS, 0.5, 0.0, 1.0, "sprint"),
        lambda: plan.lin_piece(start, (1.0, 0.0, 0.0), 1.0),
        lambda: plan.rest(start, -1.0),
        lambda: plan.Start(t=math.inf),
    ]:
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(TypeError, match="Start"):
        plan.line((5.0, 1.0, ORIGIN), (9.0, 2.0, 3.0), 1.0)


def test_orbit_entry_is_smooth():
    """Entering an orbit: the velocity is continuous where the run meets the circle."""
    rng = random.Random(43)
    centre = (1.0, -2.0, 0.5)
    outside = inside = 0
    for _ in range(300):
        start = plan.Start(t=100.0, x=(rng.uniform(-30, 30), rng.uniform(-30, 30), rng.uniform(-3, 3)))
        r, beta, a = rng.uniform(0.5, 8.0), rng.uniform(0.05, 0.9), rng.uniform(0.2, 3.0)
        try:
            entry = plan.orbit_entry(start, centre, r, beta, a)
        except ValueError:
            continue
        inside += entry.from_inside
        outside += not entry.from_inside
        circle = plan.circ_piece(entry.end, entry.orbit)
        assert circle.phase == "orbit" and circle.dtau == INF and circle.segment.t1 == INF
        w = Worldline([*entry.segments, circle.segment], anchor_t=100.0)
        w.check_continuity(1e-8 * (1.0 + math.hypot(*start.x)), 1e-8)
        # on the circle: at the radius, at the orbit's speed, counter-clockwise in (e1, e2)
        x, v = circle.segment.at(circle.segment.t0 + 0.3)
        assert math.dist(x, centre) == pytest.approx(r, rel=1e-9)
        assert math.hypot(*v) == pytest.approx(beta, abs=1e-12)
        elements = entry.orbit
        assert elements.centre == centre and elements.r == r
        assert elements.orbit.beta == beta and elements.omega == pytest.approx(beta / r)
        normal = _cross(elements.e1, elements.e2)
        radius = tuple(p - c for p, c in zip(x, centre))
        assert sum(n * c for n, c in zip(normal, _cross(radius, v))) > 0.0
        # and leaving it brakes to rest
        leave = plan.orbit_exit(plan.Start(x=x), v, a)
        assert len(leave.pieces) == 1 and math.hypot(*leave.end.v) < 1e-12
        assert leave.dphi == pytest.approx(math.atanh(beta), rel=1e-12)
    assert outside > 50 and inside > 5


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def test_planes_and_orbit_elements():
    p = plan.plane_through(ORIGIN, (3.0, 4.0, 0.0))
    assert p.e1 == pytest.approx((0.6, 0.8, 0.0)) and p.normal == pytest.approx((0.0, 0.0, 1.0))
    assert p.e2 == pytest.approx(_cross(p.normal, p.e1))
    assert repr(p).startswith("Plane(e1=(")
    # at the centre the default direction is used
    assert plan.plane_through(ORIGIN, ORIGIN).e1 == X_AXIS
    assert plan.plane_through(ORIGIN, ORIGIN, (0.0, 1.0, 0.0)).e1 == (0.0, 1.0, 0.0)

    ring = plan.OrbitElements((1.0, 1.0, 0.0), 2.0, 0.25, a0=1.0)
    assert ring.orbit == srkit.CircularOrbit(2.0, 0.5) and ring.e1 == X_AXIS and ring.a0 == 1.0
    assert "OrbitElements(centre=(1.0, 1.0, 0.0), r=2.0, omega=0.25" in repr(ring)
    lap = plan.circ_piece(plan.Start(t=3.0), ring, ring.orbit.lap_tau)
    assert lap.segment.t1 == pytest.approx(3.0 + ring.orbit.period, rel=1e-15)
    assert lap.segment.with_t1(INF) == Segment.circular(3.0, INF, (1.0, 1.0, 0.0), 2.0, 0.25, a0=1.0)
    assert lap.end().x == pytest.approx(lap.segment.at(3.0)[0])  # once round
    with pytest.raises(ValueError):
        plan.OrbitElements(ORIGIN, 2.0, 0.5)  # at the speed of light
    with pytest.raises(ValueError):
        plan.orbit_entry(plan.Start(), ORIGIN, 2.0, 1.0, 1.0)
    with pytest.raises(TypeError):
        plan.circ_piece(plan.Start(), "ring")


# ------------------------------------------------------------------------------------------------- trip


def test_trip_attributes_and_states():
    trip = srkit.Trip(D_PROXIMA, G)
    assert (trip.d, trip.a, trip.cap, trip.stops, trip.round_trip) == (D_PROXIMA, G, INF, 0, False)
    assert (trip.n, trip.n_out) == (1, 1) and trip.x_end == D_PROXIMA
    assert repr(trip) == f"Trip(d=4.2465, a={G!r}, cap=inf, stops=0, round_trip=False)"
    assert trip == srkit.Trip(D_PROXIMA, G, INF, 0, False) and trip != srkit.Trip(D_PROXIMA, G, 2.0)
    fb = srkit.FlipBurn(D_PROXIMA, G)
    assert repr(fb) == f"FlipBurn(d=4.2465, a={G!r})"
    assert trip.tau == pytest.approx(fb.tau, rel=1e-14) and trip.t_end == pytest.approx(fb.t, rel=1e-14)
    assert trip.t_out == trip.t_end and trip.leg_tau == trip.tau
    assert trip.lag == pytest.approx(fb.lag, rel=1e-12) == pytest.approx(fb.t - fb.d, rel=1e-12)
    assert trip.u_end == pytest.approx(trip.t_end - D_PROXIMA) and trip.v_end == pytest.approx(trip.t_end + D_PROXIMA)
    assert trip.du == trip.u_end and trip.dv == trip.v_end
    leg = trip.leg
    assert (leg.a, leg.tc) == (G, 0.0) and leg.phc == pytest.approx(fb.phi_peak) and leg.ta == pytest.approx(fb.tau / 2)
    assert repr(leg).startswith("Leg(a=")
    assert trip.gamma_peak == pytest.approx(fb.gamma_peak) == pytest.approx(1.0 + G * D_PROXIMA / 2.0)
    assert trip.mass_ratio() == pytest.approx(fb.mass_ratio()) == pytest.approx(math.exp(2.0 * fb.phi_peak))
    assert trip.mass_ratio(0.5) == pytest.approx(math.exp(4.0 * fb.phi_peak))

    # waiting before launch and after arrival; the turnover in the middle
    before, mid, after = trip.state(-1.0), trip.state(trip.tau / 2.0), trip.state(trip.tau + 2.0)
    assert (before.t, before.x, before.phi, before.gamma) == (-1.0, 0.0, 0.0, 1.0)
    assert mid.x == pytest.approx(D_PROXIMA / 2.0) and mid.phi == pytest.approx(fb.phi_peak)
    assert mid.gamma == pytest.approx(fb.gamma_peak) and mid.dir == 1.0
    assert mid.doppler_home == pytest.approx(math.exp(-fb.phi_peak))
    assert mid.doppler_dest == pytest.approx(math.exp(fb.phi_peak))
    assert (after.x, after.phi) == (D_PROXIMA, 0.0) and after.t == pytest.approx(trip.t_end + 2.0)
    assert trip.state_clamped(trip.tau + 2.0).t == trip.t_end and trip.state_clamped(-5.0).t == 0.0
    assert repr(mid).startswith("TripState(tau=")
    for tau in [0.1, 1.0, 2.5, 3.4]:
        s = trip.state(tau)
        assert (s.u, s.v) == pytest.approx((s.t - s.x, s.t + s.x), rel=1e-12)
        # the null coordinates invert in closed form
        assert trip.tau_of_u(s.u) == pytest.approx(tau, rel=1e-12)
        assert trip.tau_of_v(s.v) == pytest.approx(tau, rel=1e-12)


def test_a_trip_matches_its_worldline():
    """The scalar ride and the 3+1 worldline built from the plan are the same motion."""
    trip = srkit.Trip(D_PROXIMA, G)
    out = plan.line(plan.Start(), (D_PROXIMA, 0.0, 0.0), G)
    ship = Worldline([Segment.rest(-1.0, 0.0, ORIGIN), *out.segments, plan.rest(out.end, INF).segment], anchor_t=0.0)
    for tau in [0.0, 0.3, 1.2, trip.tau / 2.0, 3.0, trip.tau]:
        s, state = trip.state(tau), ship.at_tau(tau)
        assert state.t == pytest.approx(s.t, rel=1e-12, abs=1e-12)
        assert state.x[0] == pytest.approx(s.x, rel=1e-12, abs=1e-12)
        assert math.hypot(*state.v) == pytest.approx(math.tanh(s.phi), abs=1e-12)


def test_trip_light_cones_lie_on_the_cone():
    rng = random.Random(44)
    for trip in [srkit.Trip(D_PROXIMA, G), srkit.Trip(30.0, 0.7, cap=3.0, stops=2, round_trip=True)]:
        for _ in range(100):
            t, along, across = rng.uniform(-5.0, trip.t_end + 5.0), rng.uniform(-5.0, trip.d + 5.0), rng.uniform(0.0, 6.0)
            for send in (False, True):
                tau = trip.light_cone_tau(t, along, across, send)
                s = trip.state_clamped(tau) if 0.0 <= tau <= trip.tau else trip.state(tau)
                separation = math.hypot(along - s.x, across)
                gap = (t - s.t) - separation if send else (s.t - t) - separation
                assert abs(gap) <= 1e-9 * (1.0 + abs(t) + trip.t_end), (trip, t, along, across, send)
        # on the line of flight the answer is the closed-form null inverse
        assert trip.light_cone_tau(1.0, 0.0) == pytest.approx(trip.tau_of_u(1.0), rel=1e-13)
        assert trip.light_cone_tau(1.0, 0.0, send=True) == pytest.approx(trip.tau_of_v(1.0), rel=1e-13)
    with pytest.raises(ValueError):
        srkit.Trip(1.0, 1.0).light_cone_tau(0.0, 0.0, -1.0)


def test_capped_trips_stopovers_and_returns():
    there = srkit.Trip(100.0, 1.0, cap=5.0)
    assert there.gamma_peak == pytest.approx(5.0) and there.leg.tc > 0.0
    assert there.state(there.leg.ta + 0.5 * there.leg.tc).gamma == pytest.approx(5.0)
    hops = srkit.Trip(100.0, 1.0, stops=3)
    assert (hops.n, hops.n_out, hops.stops) == (4, 4, 3)
    assert hops.tau == pytest.approx(4.0 * srkit.FlipBurn(25.0, 1.0).tau)
    assert hops.state(hops.leg_tau).x == pytest.approx(25.0) and hops.state(hops.leg_tau).phi == pytest.approx(0.0, abs=1e-12)
    back = srkit.Trip(10.0, 1.0, round_trip=True)
    assert (back.n, back.n_out, back.x_end) == (2, 1, 0.0) and back.t_end == pytest.approx(2.0 * back.t_out)
    homeward = back.state(1.5 * back.leg_tau)
    assert homeward.dir == -1.0 and homeward.x == pytest.approx(5.0) and homeward.doppler_home > 1.0
    # the travelling twin is younger
    assert back.tau < back.t_end and back.state(back.tau).x == 0.0
    assert back.mass_ratio() == pytest.approx(srkit.FlipBurn(10.0, 1.0).mass_ratio() ** 2)

    for bad in [
        lambda: srkit.Trip(-1.0, 1.0),
        lambda: srkit.Trip(1.0, 0.0),
        lambda: srkit.Trip(1.0, 1.0, cap=0.5),
        lambda: srkit.Trip(math.nan, 1.0),
        lambda: srkit.FlipBurn(1.0, -1.0),
        lambda: srkit.FlipBurn(math.inf, 1.0),
        lambda: srkit.Trip(1.0, 1.0).mass_ratio(0.0),
        lambda: srkit.FlipBurn(1.0, 1.0).mass_ratio(1.5),
        lambda: srkit.Trip(1.0, 1.0, stops=-1),
    ]:
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(OverflowError):
        srkit.Trip(1.0, 1.0, stops=2**32 - 1)  # one more leg than a 32-bit count holds
    with pytest.raises(OverflowError):
        srkit.Trip(1.0, 1.0, stops=2**40)
    with pytest.raises(TypeError):
        srkit.Trip(1.0, 1.0, stops=1.5)
    with pytest.raises(OverflowError):
        srkit.FlipBurn(1e300, 1e300)  # a Lorentz factor beyond a float
    assert srkit.FlipBurn(0.0, 1.0).tau == 0.0


def test_reach():
    # from rest, straight ahead: a flip-and-burn that uses the whole budget
    budget = 3.0
    along, across = srkit.reach(0.0, G, budget, 0.0)
    assert across == pytest.approx(0.0, abs=1e-12)
    assert along == pytest.approx(srkit.flip_burn_distance(budget, G), rel=1e-9)
    assert srkit.flip_burn_tau(along, G) == pytest.approx(budget, rel=1e-9)
    # from rest the reach is the same in every direction
    for heading in [0.3, 1.5, math.pi]:
        a, b = srkit.reach(0.0, G, budget, heading)
        assert math.hypot(a, b) == pytest.approx(along, rel=1e-9)
        assert math.atan2(b, a) == pytest.approx(heading, rel=1e-9)
    # already moving: the budget must at least cancel the present speed
    assert srkit.reach(2.0, 1.0, 1.0, 0.0) is None
    stop = srkit.reach(2.0, 1.0, 2.0, 0.0)
    assert stop is not None and stop[0] == pytest.approx(math.cosh(2.0) - 1.0, rel=1e-9)  # brake at once
    with pytest.raises(ValueError):
        srkit.reach(-1.0, 1.0, 1.0, 0.0)
    with pytest.raises(ValueError):
        srkit.reach(0.0, 1.0, math.inf, 0.0)


# ------------------------------------------------------------------------------------------------ orbit


def test_circular_orbits():
    ring = srkit.CircularOrbit(100.0, 0.6)
    assert ring == srkit.CircularOrbit.from_radius_speed(100.0, 0.6)
    assert repr(ring) == "CircularOrbit(r=100.0, beta=0.6)"
    assert (ring.r, ring.beta) == (100.0, 0.6) and ring.gamma == pytest.approx(1.25)
    assert ring.phi == pytest.approx(math.atanh(0.6)) and ring.omega == pytest.approx(0.006)
    assert ring.period == pytest.approx(math.tau * 100.0 / 0.6) and ring.lap_tau == pytest.approx(ring.period / 1.25)
    assert ring.acc == pytest.approx(1.25**2 * 0.36 / 100.0)
    held = srkit.CircularOrbit.from_radius_acceleration(100.0, ring.acc)
    assert held.beta == pytest.approx(0.6) and held.gamma == pytest.approx(1.25) and held.acc == ring.acc
    assert ring.sagnac() == pytest.approx(math.tau * 100.0 * 0.6 * 1.25)
    assert ring.thomas_precession_per_lap() == pytest.approx(math.tau * 0.25)  # 2 pi (gamma - 1)
    assert ring.doppler_range() == pytest.approx((0.5, 2.0))

    state = ring.state(0.5, 10.0)
    assert state.t == pytest.approx(12.5) and state.angle == pytest.approx(0.5 + 0.006 * 12.5)
    assert math.hypot(*state.x) == pytest.approx(100.0) and math.hypot(*state.v) == pytest.approx(0.6)
    assert state.v == pytest.approx(scaled(state.heading, 0.6))
    assert repr(state).startswith("OrbitState(t=")
    seg = ring.segment(0.0, INF, ORIGIN, a0=0.5)
    assert seg == Segment.circular(0.0, INF, ORIGIN, 100.0, 0.006, a0=0.5)
    assert seg.at(12.5)[0] == pytest.approx(state.x) and seg.dtau(12.5) == pytest.approx(10.0)
    # the proper time at which the pilot sees a flash from the centre, sent at t = 0: after 100 of light travel
    assert ring.light_cone_tau(0.0, 0.0, ORIGIN) == pytest.approx(100.0 / 1.25)
    assert ring.light_cone_tau(0.0, 0.0, ORIGIN, send=True) == pytest.approx(-100.0 / 1.25)
    for bad in [
        lambda: srkit.CircularOrbit(100.0, 1.0),
        lambda: srkit.CircularOrbit(100.0, 0.0),
        lambda: srkit.CircularOrbit(0.0, 0.5),
        lambda: srkit.CircularOrbit.from_radius_acceleration(1.0, -1.0),
        lambda: ring.segment(0.0, INF, ORIGIN, X_AXIS, X_AXIS),
    ]:
        with pytest.raises(ValueError):
            bad()


# ----------------------------------------------------------------------------------- function modules


def test_rocket_equation():
    assert rocket.mass_ratio(2.0) == pytest.approx(math.e**2)
    assert rocket.mass_ratio(2.0, 0.5) == pytest.approx(math.e**4)
    assert rocket.rapidity_of_mass_ratio(rocket.mass_ratio(1.7, 0.3), 0.3) == pytest.approx(1.7)
    assert rocket.propellant_per_dry_mass(1e-12) == pytest.approx(1e-12, rel=1e-9)  # no cancellation
    assert rocket.energy_per_dry_rest_energy(1.0, 1.0, 0.5) == pytest.approx(2.0 * math.expm1(1.0))
    assert rocket.sail_rapidity(0.5 * math.expm1(1.3)) == pytest.approx(1.3)
    assert rocket.kinetic_per_rest_energy(1e-8) == pytest.approx(0.5e-16, rel=1e-9)
    assert rocket.kinetic_per_rest_energy(2.0) == pytest.approx(math.cosh(2.0) - 1.0)
    assert rocket.kinetic_per_rest_energy(2.0) == kinematics.kinetic_per_rest_energy(2.0)
    for bad in [
        lambda: rocket.mass_ratio(1.0, 0.0),
        lambda: rocket.mass_ratio(1.0, 1.1),
        lambda: rocket.rapidity_of_mass_ratio(0.0),
        lambda: rocket.energy_per_dry_rest_energy(1.0, 1.0, 0.0),
        lambda: rocket.sail_rapidity(-0.1),
    ]:
        with pytest.raises(ValueError):
            bad()


def test_kepler():
    for e in [0.0, 0.3, 0.9]:
        for m in [0.0, 0.7, 3.0, 5.5]:
            big_e = kepler.eccentric_anomaly(m, e)
            assert big_e - e * math.sin(big_e) == pytest.approx(m, abs=1e-13)
            x, y = kepler.position_in_plane(big_e, e)
            nu = kepler.true_anomaly(big_e, e)
            radius = 1.0 - e * math.cos(big_e)
            assert math.hypot(x, y) == pytest.approx(radius)
            assert (x, y) == pytest.approx((radius * math.cos(nu), radius * math.sin(nu)), abs=1e-12)
    with pytest.raises(ValueError):
        kepler.eccentric_anomaly(1.0, 1.0)
    with pytest.raises(ValueError):
        kepler.true_anomaly(1.0, -0.1)


def test_kinematics_identities():
    rng = random.Random(45)
    for _ in range(200):
        beta = scaled(unit(rng), rng.uniform(0.0, 0.99))
        u = scaled(unit(rng), rng.uniform(0.0, 0.99))
        dt, dx = rng.uniform(-10, 10), point(rng, 10.0)
        g = kinematics.gamma(beta)
        speed = math.hypot(*beta)
        assert g == pytest.approx(kinematics.gamma_of_speed(speed), rel=1e-14)
        phi = kinematics.rapidity_of_speed(speed)
        assert kinematics.speed_of_rapidity(phi) == pytest.approx(speed, abs=1e-15)
        assert kinematics.speed_shortfall(phi) == pytest.approx(1.0 - speed, rel=1e-12, abs=1e-15)
        assert kinematics.kinetic_per_rest_energy(phi) == pytest.approx(g - 1.0, rel=1e-9, abs=1e-14)
        assert kinematics.doppler_of_rapidity(phi) == pytest.approx(g * (1.0 + speed), rel=1e-13)

        # a boost keeps the interval, and the boost back undoes it
        t2, x2 = kinematics.boost(beta, dt, dx)
        assert t2 * t2 - sum(c * c for c in x2) == pytest.approx(dt * dt - sum(c * c for c in dx), rel=1e-9, abs=1e-9)
        t3, x3 = kinematics.boost(scaled(beta, -1.0), t2, x2)
        assert (t3, *x3) == pytest.approx((dt, *dx), rel=1e-10, abs=1e-10)
        # the same boost by its rapidity vector
        rapidity = scaled(beta, phi / speed) if speed > 0.0 else ORIGIN
        t4, x4 = kinematics.boost_rapidity(rapidity, dt, dx)
        assert (t4, *x4) == pytest.approx((t2, *x2), rel=1e-10, abs=1e-10)

        # composing velocities never reaches c, and the relative velocity undoes it
        w = kinematics.compose_velocities(beta, u)
        assert math.hypot(*w) < 1.0
        assert kinematics.gamma(w) == pytest.approx(g * kinematics.gamma(u) * (1.0 + sum(a * b for a, b in zip(beta, u))), rel=1e-9)
        assert kinematics.relative_velocity(beta, w) == pytest.approx(u, abs=1e-9)

        # Doppler factor and aberration of a photon, against the boost of its four-momentum
        n = unit(rng)
        k_t, k_x = kinematics.boost(beta, 1.0, n)
        assert kinematics.doppler_factor(ORIGIN, beta, n) == pytest.approx(k_t, rel=1e-12)
        assert kinematics.aberrate(n, beta) == pytest.approx(scaled(k_x, 1.0 / k_t), abs=1e-12)
        assert kinematics.doppler_factor(beta, beta, n) == pytest.approx(1.0, rel=1e-12)

    assert kinematics.aberrate_angle(math.pi / 2.0, 1.0) == pytest.approx(2.0 * math.atan(math.exp(-1.0)))
    assert kinematics.null_coordinates(5.0, 2.0) == (3.0, 7.0)
    assert kinematics.from_null_coordinates(3.0, 7.0) == (5.0, 2.0)
    a, b = (0.0, ORIGIN), (5.0, (3.0, 0.0, 0.0))
    assert kinematics.interval2(*a, *b) == 16.0 and kinematics.separation(*a, *b) == "timelike"
    assert kinematics.separation(0.0, ORIGIN, 3.0, (3.0, 0.0, 0.0)) == "null"
    assert kinematics.separation(0.0, ORIGIN, 2.0, (3.0, 0.0, 0.0)) == "spacelike"
    assert kinematics.can_reach(*a, *b) and not kinematics.can_reach(*b, *a)
    for bad in [
        lambda: kinematics.gamma((1.0, 0.0, 0.0)),
        lambda: kinematics.gamma_of_speed(1.0),
        lambda: kinematics.gamma_of_speed(-0.1),
        lambda: kinematics.boost((0.0, 1.5, 0.0), 1.0, ORIGIN),
        lambda: kinematics.compose_velocities((1.0, 0.0, 0.0), ORIGIN),
        lambda: kinematics.aberrate(ORIGIN, (0.5, 0.0, 0.0)),
        lambda: kinematics.doppler_factor(ORIGIN, ORIGIN, ORIGIN),
    ]:
        with pytest.raises(ValueError):
            bad()


def test_optics():
    # the same functions under both names
    for name in ["aberrate", "aberrate_angle", "boost_rapidity", "doppler_factor", "doppler_of_rapidity"]:
        assert getattr(optics, name) is getattr(kinematics, name)

    phi = 1.2
    ahead = optics.observe(phi, 1.0, 0.0, 0.0)
    assert ahead.doppler == pytest.approx(math.exp(phi)) and (ahead.theta, ahead.theta_seen) == (0.0, 0.0)
    assert (ahead.forward, ahead.right, ahead.up) == (1.0, 0.0, 0.0)
    astern = optics.observe(phi, -1.0, 0.0, 0.0)
    assert astern.doppler == pytest.approx(math.exp(-phi)) and astern.theta == pytest.approx(math.pi)
    abeam = optics.observe(phi, 0.0, 3.0, 4.0)
    assert abeam.doppler == pytest.approx(math.cosh(phi))  # gamma (1 + beta cos 90 deg)
    assert abeam.theta == pytest.approx(math.pi / 2.0)
    assert abeam.theta_seen == pytest.approx(kinematics.aberrate_angle(math.pi / 2.0, phi))
    assert abeam.forward == pytest.approx(math.tanh(phi))  # the sky crowds forward
    assert (abeam.right, abeam.up) == pytest.approx((0.6 / math.cosh(phi), 0.8 / math.cosh(phi)))
    still = optics.observe(phi, 0.0, 3.0, 4.0, aberration=False)
    assert still.theta_seen == still.theta and still.doppler == abeam.doppler
    assert repr(ahead).startswith("Seen(doppler=")
    # against the general formula: the source's light travels along -n towards the observer
    rng = random.Random(46)
    for _ in range(50):
        n = unit(rng)
        beta = math.tanh(phi)
        seen = optics.observe(phi, n[0], n[1], n[2])
        assert seen.doppler == pytest.approx(kinematics.doppler_factor(ORIGIN, (beta, 0.0, 0.0), scaled(n, -1.0)), rel=1e-12)
        # how fast a clock at rest there seems to run is the same factor
        assert optics.playback((beta, 0.0, 0.0), n, ORIGIN) == pytest.approx(seen.doppler, rel=1e-12)
        assert optics.playback(ORIGIN, n, ORIGIN) == 1.0

    forward, up, right = optics.basis(X_AXIS)
    assert (forward, up, right) == (X_AXIS, (0.0, 0.0, 1.0), (0.0, -1.0, 0.0))
    assert optics.basis(ORIGIN)[0] == X_AXIS

    assert optics.rayleigh_resolution(2.0, 500e-9, 1e6) == pytest.approx(1.22 * 500e-9 / 2.0 * 1e6)
    assert optics.friis_received_power(10.0, 1.0, 2.0, 0.5, 4.0) == pytest.approx(10.0 * (math.pi * 2.0 / 8.0) ** 2)
    photon = srkit.units.PLANCK * srkit.units.C / 1e-6
    assert optics.photon_bit_rate(photon * 1e6, 1e-6) == pytest.approx(1e6)
    assert optics.photon_bit_rate(photon, 1e-6, 3.0) == pytest.approx(3.0)
    sun = optics.black_body_luminosity(srkit.units.R_SUN, 5772.0)
    assert sun == pytest.approx(srkit.units.L_SUN, rel=2e-3)
    for bad in [
        lambda: optics.observe(1.0, 0.0, 0.0, 0.0),
        lambda: optics.observe(math.inf, 1.0, 0.0, 0.0),
        lambda: optics.playback((1.0, 0.0, 0.0), X_AXIS, ORIGIN),
    ]:
        with pytest.raises(ValueError):
            bad()
