"""Segments: the identities that must hold whatever the numbers, and the shape of the Python API."""

from __future__ import annotations

import math
import random

import pytest

import srkit
from helpers import G, INF, ORIGIN, X_AXIS, point, scaled, segments_of_every_kind, unit
from srkit import Segment

KINDS = ["rest", "inertial", "hyperbolic", "hyperbolic", "circular", "steered", "burn"]


def cases(seed: int, rounds: int = 40):
    rng = random.Random(seed)
    for _ in range(rounds):
        yield from segments_of_every_kind(rng)


def test_the_builder_covers_every_kind_but_null():
    assert [s.kind for s in segments_of_every_kind(random.Random(0))] == KINDS


def test_proper_time_inverts_coordinate_time():
    """`t_of_tau(dtau(t)) == t` and back, for every kind with a clock."""
    rng = random.Random(11)
    for seg in cases(1):
        for _ in range(5):
            t = rng.uniform(seg.t0, seg.t1)
            s = seg.dtau(t)
            assert 0.0 <= s <= t - seg.t0 + 1e-12, (seg, t)
            assert seg.t_of_tau(s) == pytest.approx(t, rel=1e-11, abs=1e-11), (seg, t)
            assert seg.dtau(seg.t_of_tau(s)) == pytest.approx(s, rel=1e-11, abs=1e-11), (seg, t)


def test_the_state_by_proper_time_is_the_state_by_coordinate_time():
    rng = random.Random(12)
    for seg in cases(2):
        t = rng.uniform(seg.t0, seg.t1)
        x, v = seg.at(t)
        state = seg.at_tau(seg.dtau(t))
        assert state.t == pytest.approx(t, rel=1e-11, abs=1e-11)
        assert math.dist(state.x, x) <= 1e-9 * (1.0 + math.hypot(*x)), seg
        assert math.dist(state.v, v) <= 1e-9, seg
        assert state.jump is False
        assert math.hypot(*v) < 1.0
        # a clock on anything that moves runs slow: d tau / dt = 1 / gamma
        h = 1e-6
        rate = (seg.dtau(t + h) - seg.dtau(t - h)) / (2.0 * h)
        assert rate == pytest.approx(math.sqrt(1.0 - math.hypot(*v) ** 2), rel=1e-5, abs=1e-7), seg


def test_proper_time_is_monotonic_and_starts_at_zero():
    for seg in cases(3, rounds=10):
        assert seg.dtau(seg.t0) == 0.0
        ts = [seg.t0 + (seg.t1 - seg.t0) * k / 20.0 for k in range(21)]
        taus = [seg.dtau(t) for t in ts]
        assert all(b > a for a, b in zip(taus, taus[1:])), seg
        assert seg.proper_duration() == taus[-1]


def test_light_cone_roots_lie_on_the_cone():
    """`recv` gives the moment light from an event arrives: t - tE == |x(t) - xE|. `send` gives the moment
    seen from an event: tE - t == |xE - x(t)|."""
    rng = random.Random(13)
    received = sent = 0
    for seg in cases(4):
        for _ in range(4):
            te, xe = rng.uniform(seg.t0 - 60.0, seg.t1 + 60.0), point(rng)
            scale = 1.0 + abs(te) + math.hypot(*xe)
            try:
                t = seg.recv(te, xe, seg.t0, seg.t1)
            except srkit.NoIntersection:
                pass
            else:
                received += 1
                assert seg.t0 <= t <= seg.t1
                gap = (t - te) - math.dist(seg.at(t)[0], xe)
                # the light had already passed at t0 (the answer is then t0), or the root is on the cone
                assert (t == seg.t0 and gap >= 0.0) or abs(gap) <= 1e-9 * scale, (seg, te, xe, gap)
            try:
                t = seg.send(te, xe, seg.t0, seg.t1)
            except srkit.NoIntersection:
                continue
            sent += 1
            assert seg.t0 <= t <= seg.t1
            gap = (te - t) - math.dist(xe, seg.at(t)[0])
            assert abs(gap) <= 1e-9 * scale, (seg, te, xe, gap)
    # the sample must exercise both queries
    assert received > 100 and sent > 100


def test_recv_and_send_are_inverse_queries():
    """Light sent from the segment at t reaches an observer at (tE, xE); light from an event at the segment's
    own position at t arrives there at t."""
    rng = random.Random(14)
    for seg in cases(5, rounds=15):
        t = rng.uniform(seg.t0, seg.t1)
        x, _ = seg.at(t)
        xe = point(rng)
        te = t + math.dist(x, xe)  # where and when the light from (t, x) is at xe
        assert seg.send(te, xe) == pytest.approx(t, rel=1e-9, abs=1e-9), seg
        te = t - math.dist(x, xe)  # an event whose light reaches x at t
        assert seg.recv(te, xe) == pytest.approx(t, rel=1e-9, abs=1e-9), seg


def test_null_segments_and_jumps():
    n = (0.6, 0.0, 0.8)
    ray = Segment.null(2.0, 7.0, (1.0, 1.0, 1.0), n)
    assert ray.kind == "null" and ray.is_null
    assert ray.dtau(5.0) == 0.0 and ray.proper_duration() == 0.0 and ray.vmax() == 1.0
    x, v = ray.at(4.0)
    assert x == pytest.approx((2.2, 1.0, 2.6)) and v == n
    assert ray.t_of_tau(0.0) == 7.0
    after = ray.at_tau(0.0)
    assert after.jump is True and after.t == 7.0 and after.x == pytest.approx((4.0, 1.0, 5.0))
    assert after.v == (0.0, 0.0, 0.0)

    jump = Segment.jump(5.0, (1.0, 2.0, 3.0), (1.0, 2.0, 7.0))
    assert jump.t1 == 9.0 and jump.at_tau(0.0).x == (1.0, 2.0, 7.0)
    assert jump.params["p1"] == (1.0, 2.0, 7.0) and jump.params["n"] == (0.0, 0.0, 1.0)
    # light coming the other way from a point ahead meets the jump part of the way along
    assert jump.recv(5.0, (1.0, 2.0, 6.0), 5.0, 9.0) == pytest.approx(6.5)
    # light that sets off behind the jump after it has left never catches it
    with pytest.raises(srkit.NoIntersection):
        jump.recv(4.5, (1.0, 2.0, 2.0), 5.0, 9.0)
    # light that runs alongside from the start: the first moment it is there
    assert jump.recv(4.0, (1.0, 2.0, 2.0), 5.0, 9.0) == 5.0
    # a point off to the side sees the jump pass
    seen = jump.send(12.0, (4.0, 2.0, 5.0), 5.0, 9.0)
    assert (12.0 - seen) == pytest.approx(math.dist((4.0, 2.0, 5.0), jump.at(seen)[0]))


def test_hyperbolic_closed_forms():
    a = G
    seg = Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, a)
    for t in [1e-9, 0.3, 1.0, 12.0, 1e5, 1e100]:
        x, v = seg.at(t)
        gamma = math.hypot(1.0, a * t)
        assert v[0] == pytest.approx(a * t / gamma, rel=1e-15)
        assert seg.dtau(t) == pytest.approx(math.asinh(a * t) / a, rel=1e-14)
        if t < 1e50:
            # x = (gamma - 1) / a, formed without cancellation
            assert x[0] == pytest.approx((a * t) ** 2 / (gamma + 1.0) / a, rel=1e-14)
    assert seg.vmax() == 1.0
    assert seg.with_t1(2.0).vmax() == pytest.approx(math.tanh(math.asinh(2.0 * a)))
    state = seg.at_tau(3.0)
    assert state.phi == pytest.approx(3.0 * a) and state.acc == a
    assert state.t == pytest.approx(math.sinh(3.0 * a) / a)
    # the default rapidity at the start is zero
    assert seg == Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, a, 0.0)
    assert seg.params == {"p0": ORIGIN, "n": X_AXIS, "alpha": a, "phi0": 0.0}


def test_a_steered_burn_along_the_motion_is_a_hyperbola():
    rng = random.Random(15)
    for _ in range(50):
        n, p0 = unit(rng), point(rng)
        a, phi0 = rng.uniform(0.05, 3.0), rng.uniform(0.0, 3.0)
        hyp = Segment.hyperbolic(1.0, 9.0, p0, n, a, phi0)
        steer = Segment.steered(1.0, 9.0, p0, scaled(n, math.tanh(phi0)), n, a)
        proper = Segment.steered(1.0, 9.0, p0, scaled(n, math.sinh(phi0)), n, a, proper_velocity=True)
        for t in [1.0, 2.5, 9.0]:
            for other in (steer, proper):
                assert math.dist(other.at(t)[0], hyp.at(t)[0]) <= 1e-9 * (1.0 + math.hypot(*p0))
                assert math.dist(other.at(t)[1], hyp.at(t)[1]) <= 1e-9
                assert other.dtau(t) == pytest.approx(hyp.dtau(t), rel=1e-9, abs=1e-12)
    # thrust at right angles from rest: the ship still follows a hyperbola, along the thrust
    side = Segment.steered(0.0, 5.0, ORIGIN, ORIGIN, (0.0, 2.0, 0.0), 1.0)
    assert side.params["thrust"] == (0.0, 1.0, 0.0)  # the thrust direction is normalised
    assert side.at(3.0)[0] == pytest.approx((0.0, math.hypot(1.0, 3.0) - 1.0, 0.0))


def test_a_burn_obeys_the_rocket_equation():
    w, alpha0, phi0 = 0.6, 0.4, 0.2
    t_burn = w / alpha0
    seg = Segment.burn(0.0, 1.0, ORIGIN, X_AXIS, phi0, alpha0, w)
    assert seg.kind == "burn"
    assert seg.params["t_burn"] == pytest.approx(t_burn) and seg.params["alpha0"] == pytest.approx(alpha0)
    for s in [0.0, 0.1, 0.5, 1.0]:
        state = seg.at_tau(s)
        mass_left = 1.0 - s / t_burn
        # rapidity gained = w ln(m0 / m)
        assert state.phi == pytest.approx(phi0 + w * math.log(1.0 / mass_left), rel=1e-13)
        assert state.phi - phi0 == pytest.approx(
            srkit.rocket.rapidity_of_mass_ratio(1.0 / mass_left, w), rel=1e-12, abs=1e-15
        )
        # the push is constant, so the acceleration rises as the ship lightens
        assert state.acc == pytest.approx(alpha0 / mass_left, rel=1e-13)
    # braking burns thrust against n
    brake = Segment.burn(0.0, 1.0, ORIGIN, X_AXIS, 1.0, -alpha0, w)
    assert brake.at_tau(0.5).phi < 1.0 and brake.at_tau(0.5).acc < 0.0
    # the mass would be gone at proper time w / alpha0: the segment has to end first
    with pytest.raises(ValueError, match="mass runs out"):
        Segment.burn(0.0, 1e6, ORIGIN, X_AXIS, phi0, alpha0, w)
    with pytest.raises(ValueError, match="mass runs out"):
        Segment.burn(0.0, INF, ORIGIN, X_AXIS, phi0, alpha0, w)
    with pytest.raises(ValueError, match="mass runs out"):
        seg.with_t1(1e6)
    # before its start a burn continues backwards as a formula (heavier and slower, no NaN); light-cone
    # queries keep to the span where the ship really is burning
    before = seg.at(-1.0)
    assert all(math.isfinite(c) for c in before[0]) and seg.dtau(-1.0) < 0.0
    assert before[1][0] < seg.at(0.0)[1][0]
    burnt_out = seg.t_of_tau(t_burn * (1.0 - 2.0**-52))
    answered = 0
    for te, xe in [(-5.0, (3.0, 0.0, 0.0)), (0.5, (0.2, 0.1, 0.0)), (1.5, (0.9, 0.0, 0.3)), (50.0, (1.0, 0.0, 0.0))]:
        for query in (seg.recv, seg.send):
            try:
                t = query(te, xe)
            except srkit.NoIntersection:
                continue
            assert 0.0 <= t <= burnt_out
            answered += 1
    assert answered >= 2


def test_circular_motion():
    r, beta = 3.0, 0.8
    seg = Segment.circular(1.0, 50.0, (1.0, 2.0, 3.0), r, beta / r, a0=0.5)
    assert seg.params["e1"] == X_AXIS and seg.params["e2"] == (0.0, 1.0, 0.0)
    for t in [1.0, 4.0, 33.0]:
        x, v = seg.at(t)
        assert math.dist(x, (1.0, 2.0, 3.0)) == pytest.approx(r)
        assert math.hypot(*v) == pytest.approx(beta)
        assert seg.dtau(t) == pytest.approx((t - 1.0) * 0.6)
    assert seg.vmax() == pytest.approx(beta)
    state = seg.at_tau(2.0)
    assert state.phi == pytest.approx(math.atanh(beta))
    assert state.acc == pytest.approx(beta**2 / (1.0 - beta**2) / r)  # gamma^2 beta^2 / r
    # any orthonormal plane will do
    tilted = Segment.circular(0.0, 1.0, ORIGIN, r, beta / r, (0.0, 0.0, 1.0), (0.0, -1.0, 0.0))
    assert tilted.at(0.0)[0] == pytest.approx((0.0, 0.0, r))


def test_horizons():
    a = 0.5
    ship = Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, a)
    horizon = ship.rindler_horizon()
    assert horizon.distance == 1.0 / a and horizon.n == X_AXIS
    assert horizon.pivot_t == 0.0 and horizon.pivot_x == (-1.0 / a, 0.0, 0.0)
    # a point 1.5 behind the start: its signals stop arriving after t = 0.5
    behind = (-1.5, 0.0, 0.0)
    assert horizon.last_signal_time(behind) == pytest.approx(0.5)
    assert horizon.can_signal(0.25, behind) and not horizon.can_signal(0.75, behind)
    assert horizon.margin(0.25, behind) == pytest.approx(0.25)
    assert ship.horizon_margin(0.25, behind) == pytest.approx(0.25)
    assert ship.horizon_margin(0.75, behind) == pytest.approx(-0.25)
    assert ship.horizon_margin(0.75, behind, past=True) == pytest.approx(1.25)
    assert horizon.margin_to_be_seen(0.75, behind) == pytest.approx(1.25)
    # inside the horizon the light arrives; beyond it, never
    t = ship.recv(0.25, behind)
    assert (t - 0.25) == pytest.approx(math.dist(ship.at(t)[0], behind), rel=1e-12)
    with pytest.raises(srkit.NoIntersection):
        ship.recv(0.75, behind)
    # a steered segment has the same horizon; the other kinds have none
    steered = Segment.steered(0.0, INF, ORIGIN, ORIGIN, X_AXIS, a)
    assert steered.horizon_margin(0.25, behind) == pytest.approx(0.25)
    assert steered.rindler_horizon() is None
    rest = Segment.rest(0.0, INF, ORIGIN)
    assert rest.horizon_margin(0.0, behind) is None and rest.rindler_horizon() is None
    assert "RindlerHorizon(pivot_t=0.0" in repr(horizon)


def test_no_intersection_outside_the_range():
    rest = Segment.rest(0.0, 10.0, ORIGIN)
    assert rest.recv(1.0, (3.0, 0.0, 0.0)) == 4.0
    assert rest.recv(1.0, (3.0, 0.0, 0.0), 0.0, 10.0) == 4.0
    with pytest.raises(srkit.NoIntersection):
        rest.recv(1.0, (3.0, 0.0, 0.0), 0.0, 3.5)
    # the light had already passed at lo: the answer is lo
    assert rest.recv(1.0, (3.0, 0.0, 0.0), 6.0, 10.0) == 6.0
    assert rest.send(9.0, (3.0, 0.0, 0.0)) == 6.0
    with pytest.raises(srkit.NoIntersection):
        rest.send(9.0, (3.0, 0.0, 0.0), 7.0, 10.0)
    with pytest.raises(srkit.NoIntersection):
        rest.send(9.0, (3.0, 0.0, 0.0), 0.0, 5.0)
    assert issubclass(srkit.NoIntersection, Exception)
    assert not issubclass(srkit.NoIntersection, ValueError)


class Triple:
    """A minimal sequence of three numbers that is neither a tuple nor a list."""

    def __init__(self, *values):
        self.values = values

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index):
        return self.values[index]


def test_vectors_are_numpy_arrays_too():
    np = pytest.importorskip("numpy")
    expected = Segment.inertial(0.0, 1.0, (1.0, 2.0, 3.0), (0.1, 0.2, 0.3))
    forms = [
        lambda v: np.array(v),
        lambda v: np.array([v, v])[1],  # a row of a larger array
        lambda v: [np.float64(c) for c in v],
    ]
    for form in forms:
        assert Segment.inertial(0.0, 1.0, form((1.0, 2.0, 3.0)), form((0.1, 0.2, 0.3))) == expected
    for p in (np.array([1, 2, 3]), np.array([1, 2, 3], dtype=np.float32)):
        assert Segment.inertial(0.0, 1.0, p, (0.1, 0.2, 0.3)) == expected
    assert expected.recv(0.0, np.array([5.0, 0.0, 0.0])) == expected.recv(0.0, (5.0, 0.0, 0.0))
    x, v = expected.at(np.float64(0.5))
    assert type(x) is tuple and all(type(c) is float for c in x + v)
    with pytest.raises(ValueError, match="three components"):
        Segment.rest(0.0, 1.0, np.zeros(4))
    with pytest.raises(ValueError, match="three components"):
        Segment.rest(0.0, 1.0, np.zeros((3, 3)).ravel())


def test_vectors_are_any_sequence_of_three_numbers():
    expected = Segment.inertial(0.0, 1.0, (1.0, 2.0, 3.0), (0.1, 0.2, 0.3))
    for form in (list, lambda v: Triple(*v), tuple):
        assert Segment.inertial(0.0, 1.0, form((1.0, 2.0, 3.0)), form((0.1, 0.2, 0.3))) == expected
    # whole numbers in any form
    for p in ([1, 2, 3], range(1, 4), (True, 2, 3.0)):
        assert Segment.inertial(0.0, 1.0, p, (0.1, 0.2, 0.3)) == expected
    # ints are numbers too, and what comes back is always a tuple of floats
    x, v = Segment.rest(0, 1, [1, 2, 3]).at(0)
    assert x == (1.0, 2.0, 3.0) and v == (0.0, 0.0, 0.0)
    assert all(type(c) is float for c in x + v)

    with pytest.raises(ValueError, match="three components"):
        Segment.rest(0.0, 1.0, (1.0, 2.0))
    with pytest.raises(ValueError, match="three components"):
        Segment.rest(0.0, 1.0, [1.0, 2.0, 3.0, 4.0])
    with pytest.raises(TypeError):
        Segment.rest(0.0, 1.0, 5.0)
    with pytest.raises(TypeError):
        Segment.rest(0.0, 1.0, "abc")
    with pytest.raises(TypeError):
        Segment.rest(0.0, 1.0, (1.0, "b", 3.0))
    with pytest.raises(TypeError):
        Segment.rest("0", 1.0, ORIGIN)


def test_constructors_refuse_what_is_not_physics():
    bad = [
        lambda: Segment.rest(0.0, math.nan, ORIGIN),
        lambda: Segment.rest(2.0, 1.0, ORIGIN),
        lambda: Segment.rest(0.0, 1.0, (math.inf, 0.0, 0.0)),
        lambda: Segment.inertial(0.0, 1.0, ORIGIN, (1.0, 0.0, 0.0)),  # at the speed of light
        lambda: Segment.inertial(0.0, 1.0, ORIGIN, (0.8, 0.8, 0.0)),  # faster than light
        lambda: Segment.inertial(0.0, 1.0, ORIGIN, (math.nan, 0.0, 0.0)),
        lambda: Segment.hyperbolic(0.0, 1.0, ORIGIN, (2.0, 0.0, 0.0), 1.0),  # n is not a unit vector
        lambda: Segment.hyperbolic(0.0, 1.0, ORIGIN, ORIGIN, 1.0),
        lambda: Segment.hyperbolic(0.0, 1.0, ORIGIN, X_AXIS, 0.0),
        lambda: Segment.hyperbolic(0.0, 1.0, ORIGIN, X_AXIS, math.inf),
        lambda: Segment.hyperbolic(0.0, 1.0, ORIGIN, X_AXIS, 1.0, math.nan),
        lambda: Segment.steered(0.0, 1.0, ORIGIN, (1.0, 0.0, 0.0), X_AXIS, 1.0),
        lambda: Segment.steered(0.0, 1.0, ORIGIN, ORIGIN, ORIGIN, 1.0),  # no thrust direction
        lambda: Segment.steered(0.0, 1.0, ORIGIN, ORIGIN, X_AXIS, -1.0),
        lambda: Segment.burn(0.0, 0.1, ORIGIN, X_AXIS, 0.0, 0.0, 1.0),
        lambda: Segment.burn(0.0, 0.1, ORIGIN, X_AXIS, 0.0, 1.0, 1.5),  # exhaust faster than light
        lambda: Segment.burn(0.0, 0.1, ORIGIN, X_AXIS, 0.0, 1.0, 0.0),
        lambda: Segment.circular(0.0, 1.0, ORIGIN, 2.0, 0.5),  # r omega = 1
        lambda: Segment.circular(0.0, 1.0, ORIGIN, -2.0, 0.1),
        lambda: Segment.circular(0.0, 1.0, ORIGIN, 2.0, 0.1, X_AXIS, X_AXIS),  # not a plane
        lambda: Segment.circular(0.0, 1.0, ORIGIN, 2.0, 0.1, (3.0, 0.0, 0.0)),
        lambda: Segment.null(0.0, 1.0, ORIGIN, (1.0, 1.0, 0.0)),
        lambda: Segment.null(0.0, INF, ORIGIN, X_AXIS),
        lambda: Segment.jump(INF, ORIGIN, X_AXIS),
        lambda: Segment.rest(0.0, 1.0, ORIGIN).with_t1(-1.0),
        # a moving segment is pinned at its start, which must then be a finite time
        lambda: Segment.inertial(-INF, 1.0, ORIGIN, (0.5, 0.0, 0.0)),
        lambda: Segment.hyperbolic(-INF, 1.0, ORIGIN, X_AXIS, 1.0),
        lambda: Segment.circular(-INF, 1.0, ORIGIN, 2.0, 0.1),
        lambda: Segment.steered(-INF, 1.0, ORIGIN, ORIGIN, X_AXIS, 1.0),
    ]
    for make in bad:
        with pytest.raises(ValueError):
            make()
    # a unit vector to within rounding is accepted and passed through bit for bit
    n = (1.0 / math.sqrt(3.0),) * 3
    assert Segment.hyperbolic(0.0, 1.0, ORIGIN, n, 1.0).params["n"] == n
    # a rest may reach back for ever, and any segment on for ever
    assert Segment.rest(-INF, INF, ORIGIN).t0 == -INF
    assert Segment.inertial(0.0, INF, ORIGIN, (0.5, 0.0, 0.0)).t1 == INF
    assert Segment.rest(0.0, 0.0, ORIGIN).t1 == 0.0


def test_repr_reads_back():
    namespace = {"Segment": Segment, "inf": INF}
    rng = random.Random(16)
    for seg in [*segments_of_every_kind(rng), Segment.rest(0.0, INF, ORIGIN)]:
        again = eval(repr(seg), namespace)
        assert again.kind == seg.kind and again.t0 == seg.t0 and again.t1 == seg.t1
        if seg.kind in ("rest", "inertial", "hyperbolic", "circular"):
            assert again == seg, repr(seg)
        for t in [seg.t0, 0.5 * (seg.t0 + min(seg.t1, seg.t0 + 10.0))]:
            assert math.dist(again.at(t)[0], seg.at(t)[0]) <= 1e-9
    ray = Segment.null(0.0, 2.0, ORIGIN, (0.0, 1.0, 0.0))
    assert repr(ray) == "Segment.null(t0=0.0, t1=2.0, p0=(0.0, 0.0, 0.0), n=(0.0, 1.0, 0.0))"
    assert eval(repr(ray), namespace) == ray
    assert repr(Segment.rest(0.0, INF, ORIGIN)) == "Segment.rest(t0=0.0, t1=inf, p=(0.0, 0.0, 0.0))"
    state = Segment.rest(1.0, 2.0, ORIGIN).at_tau(0.5)
    assert repr(state) == (
        "TauState(t=1.5, x=(0.0, 0.0, 0.0), v=(0.0, 0.0, 0.0), phi=0.0, acc=0.0, jump=False)"
    )


def test_segments_are_values():
    a = Segment.inertial(0.0, 1.0, ORIGIN, (0.5, 0.0, 0.0))
    b = Segment.inertial(0.0, 1.0, ORIGIN, (0.5, 0.0, 0.0))
    assert a == b and a is not b and not (a != b)
    assert a != a.with_t1(2.0) and a.with_t1(2.0).t1 == 2.0 and a.t1 == 1.0
    assert a != Segment.rest(0.0, 1.0, ORIGIN) and a != "a segment"
    with pytest.raises(TypeError):
        hash(a)
    with pytest.raises(TypeError):
        a < b
    # there is nothing to assign to: a segment cannot be changed
    with pytest.raises(AttributeError):
        a.t1 = 5.0
    with pytest.raises(TypeError):
        Segment()
    # params is a fresh dict each time
    a.params["v"] = None
    assert a.params == {"p0": ORIGIN, "v": (0.5, 0.0, 0.0)}
