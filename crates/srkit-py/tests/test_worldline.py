"""Worldlines: proper time along the chain, light-cone queries, and the vectorised methods."""

from __future__ import annotations

import array
import math
import random

import pytest

import srkit
from helpers import D_PROXIMA, G, INF, ORIGIN, X_AXIS, point, proxima_ship
from srkit import Segment, Worldline, plan

HOME = ORIGIN
PROXIMA = (D_PROXIMA, 0.0, 0.0)


def tour() -> Worldline:
    """A worldline with every kind of segment: wait, boost at an angle, coast, burn, coast, jump, orbit."""
    start = plan.Start(t=10.0, tau=3.0, x=(1.0, -2.0, 0.5))
    out = plan.line(start, (6.0, 2.0, 0.5), a=1.5, phi_cap=math.atanh(0.8))
    entry = plan.orbit_entry(out.end, (20.0, 5.0, 0.5), 3.0, 0.6, a=1.2)
    circle = plan.circ_piece(entry.end, entry.orbit, 40.0)
    leave = plan.orbit_exit(circle.end(), circle.end().v, 2.0)
    hop = plan.null_piece(leave.end, (30.0, 30.0, 3.0))
    segments = [
        Segment.rest(0.0, 10.0, start.x),
        *out.segments,
        *entry.segments,
        circle.segment,
        *leave.segments,
        hop.segment,
        plan.rest(hop.end(), INF).segment,
    ]
    return Worldline(segments, anchor_t=10.0, anchor_tau=3.0)


def sample_times(w: Worldline, rng: random.Random, n: int = 200) -> list[float]:
    last = w[-1].t0
    return sorted(rng.uniform(w[0].t0 - 5.0, last + 20.0) for _ in range(n))


def test_a_worldline_is_a_sequence_of_segments():
    w = tour()
    assert len(w) == len(w.segments) == len(list(w)) == len(w.tau0s)
    assert all(isinstance(seg, Segment) for seg in w)
    assert w[0] == w.segments[0] and w[-1] == w.segments[-1] and w[len(w) - 1] == w[-1]
    assert list(w) == list(w.segments)
    with pytest.raises(IndexError):
        w[len(w)]
    with pytest.raises(IndexError):
        w[-len(w) - 1]
    kinds = [seg.kind for seg in w]
    assert kinds[0] == "rest" and kinds[-1] == "rest"
    assert {"hyperbolic", "inertial", "circular", "null"} <= set(kinds)
    assert w.has_null and not proxima_ship().has_null
    assert repr(proxima_ship()) == (
        "Worldline(4 segments: rest, hyperbolic, hyperbolic, rest; t0=-1.0, tau0=-1.0)"
    )
    assert isinstance(w, Worldline) and type(w) is Worldline
    w.check_continuity(1e-9, 1e-9)
    assert w.is_continuous()
    assert w.vmax() == 1.0  # the jump moves at c
    assert proxima_ship().vmax() == pytest.approx(math.tanh(srkit.FlipBurn(D_PROXIMA, G).phi_peak))


def test_proper_time_runs_along_the_chain():
    w = tour()
    assert w.tau(10.0) == 3.0
    assert w.tau0s[0] == pytest.approx(-7.0)  # at rest for ten years before the anchor
    rng = random.Random(21)
    previous = -INF
    for t in sample_times(w, rng):
        tau = w.tau(t)
        assert tau >= previous
        previous = tau
        i = w.find(t)
        assert w[i].t0 <= t or i == 0
        assert tau == pytest.approx(w.tau0s[i] + w[i].dtau(t), rel=1e-15, abs=1e-15)
        assert w.at(t) == w[i].at(t)
        if w[i].kind != "null":
            # and back: the same moment by proper time
            assert w.t_of_tau(tau) == pytest.approx(t, rel=1e-10, abs=1e-10)
            state = w.at_tau(tau)
            assert state.t == pytest.approx(t, rel=1e-10, abs=1e-10)
            assert math.dist(state.x, w.at(t)[0]) <= 1e-8
            assert w.find_tau(tau) in (i, i + 1, i - 1)
    # no proper time passes on the jump
    j = [seg.kind for seg in w].index("null")
    assert w.tau(w[j].t0) == pytest.approx(w.tau(w[j + 1].t0), abs=1e-12)


def test_the_anchor_sets_the_clock():
    segments = proxima_ship().segments
    ship = Worldline(segments, anchor_t=0.0)
    late = Worldline(segments, anchor_t=2.0, anchor_tau=100.0)
    assert ship.tau(0.0) == 0.0 and late.tau(2.0) == 100.0
    # by default the clock reads zero at the start of the first segment
    assert Worldline(segments).tau(-1.0) == 0.0 and Worldline(segments).tau(0.0) == 1.0
    for t in [-3.0, 0.5, 2.0, 4.0, 9.0]:
        assert late.tau(t) - ship.tau(t) == pytest.approx(100.0 - ship.tau(2.0), rel=1e-13)
    # the ride takes the flip-and-burn's proper time, and the traveller then ages at the home rate
    fb = srkit.FlipBurn(D_PROXIMA, G)
    assert ship.tau(fb.t) == pytest.approx(fb.tau, rel=1e-13)
    assert ship.tau(fb.t + 5.0) == pytest.approx(fb.tau + 5.0, rel=1e-13)
    assert ship.at(fb.t + 5.0) == (PROXIMA, ORIGIN)


def test_worldline_light_cones_agree_with_the_segments():
    """`Worldline.recv` and `send` find the same roots as the segment-level calls on the segment in force."""
    w = tour()
    rng = random.Random(22)
    checked = 0
    for _ in range(400):
        te, xe = rng.uniform(-20.0, 120.0), point(rng, 40.0)
        scale = 1.0 + abs(te) + math.hypot(*xe)

        t = w.recv(te, xe)
        assert math.isfinite(t) and t >= te
        x, _ = w.at(t)
        assert abs((t - te) - math.dist(x, xe)) <= 1e-9 * scale
        i = w.find(t)
        lo = -INF if i == 0 else w[i].t0
        hi = INF if i == len(w) - 1 else w[i + 1].t0
        try:
            assert w[i].recv(te, xe, lo, hi) == pytest.approx(t, rel=1e-12, abs=1e-12)
            checked += 1
        except srkit.NoIntersection:
            # the root sits on a join and belongs to the neighbour by a rounding
            assert min(abs(t - lo), abs(t - hi)) <= 1e-9 * scale

        t = w.send(te, xe)
        assert math.isfinite(t) and t <= te
        x, _ = w.at(t)
        assert abs((te - t) - math.dist(xe, x)) <= 1e-9 * scale
        i = w.find(t)
        lo = -INF if i == 0 else w[i].t0
        hi = INF if i == len(w) - 1 else w[i + 1].t0
        try:
            assert w[i].send(te, xe, lo, hi) == pytest.approx(t, rel=1e-12, abs=1e-12)
            checked += 1
        except srkit.NoIntersection:
            assert min(abs(t - lo), abs(t - hi)) <= 1e-9 * scale
    assert checked > 700


def test_what_home_sees_of_the_ship():
    ship = proxima_ship()
    fb = srkit.FlipBurn(D_PROXIMA, G)
    # the arrival is seen at home one light-crossing after it happens
    assert ship.send(fb.t + D_PROXIMA, HOME) == pytest.approx(fb.t, rel=1e-13)
    # before the launch the ship is seen where it was: at home, now
    assert ship.send(-2.0, HOME) == -2.0
    t_emit, x, v = ship.seen_from(3.0, HOME)
    assert t_emit == ship.send(3.0, HOME) and (x, v) == ship.at(t_emit)
    assert 3.0 - t_emit == pytest.approx(x[0], rel=1e-13)  # light took as long as the ship was far
    # a message sent from home half a year after the launch catches the ship while it is still under way
    t_recv = ship.recv(0.5, HOME)
    assert t_recv - 0.5 == pytest.approx(ship.at(t_recv)[0][0], rel=1e-13)
    assert 0.5 < t_recv < fb.t
    # the same moments by the trip planner's null coordinates
    trip = srkit.Trip(D_PROXIMA, G)
    assert ship.tau(t_recv) == pytest.approx(trip.tau_of_u(0.5), rel=1e-12)
    assert ship.tau(t_emit) == pytest.approx(trip.tau_of_v(3.0), rel=1e-12)


def test_an_everlasting_burn_outruns_light():
    a = 0.5
    runaway = Worldline([Segment.rest(-5.0, 0.0, ORIGIN), Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, a)])
    behind = (-1.5, 0.0, 0.0)
    assert math.isfinite(runaway.recv(0.25, behind))
    assert runaway.recv(0.75, behind) == INF  # sent from beyond the horizon
    radar = runaway.radar(0.75, behind)
    assert radar.tau_recv == INF and math.isfinite(radar.tau_send)
    assert radar.time == INF and radar.distance == INF
    with pytest.raises(srkit.NoIntersection):
        runaway.light_cone_coordinates(0.75, behind)
    # nothing of the runaway is hidden from behind, though: it can always be seen
    assert math.isfinite(runaway.send(50.0, behind))
    # a ship that has always been coming in at speed cannot be seen from behind its past horizon
    incoming = Worldline([Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, a)])
    assert incoming.send(-3.0, (-1.0, 0.0, 0.0)) == -INF
    with pytest.raises(srkit.NoIntersection):
        incoming.seen_from(-3.0, (-1.0, 0.0, 0.0))


def test_radar_and_light_cone_coordinates():
    home = Worldline([Segment.rest(0.0, INF, ORIGIN)])
    radar = home.radar(10.0, (3.0, 4.0, 0.0))
    assert (radar.tau_send, radar.tau_recv) == (5.0, 15.0)
    assert radar.time == 10.0 and radar.distance == 5.0
    assert repr(radar) == "Radar(tau_send=5.0, tau_recv=15.0)"
    assert home.light_cone_coordinates(10.0, (3.0, 4.0, 0.0)) == (15.0, 5.0)
    # a moving observer: radar time on the observer's own worldline is its proper time
    ship = proxima_ship()
    for t in [0.5, 2.0, 5.0]:
        x, _ = ship.at(t)
        own = ship.radar(t, x)
        assert own.time == pytest.approx(ship.tau(t), rel=1e-12) and abs(own.distance) < 1e-12


def test_continuity_is_checked_not_assumed():
    good = [Segment.rest(0.0, 1.0, ORIGIN), Segment.hyperbolic(1.0, INF, ORIGIN, X_AXIS, 1.0)]
    Worldline(good).check_continuity()
    gap = [Segment.rest(0.0, 1.0, ORIGIN), Segment.hyperbolic(1.0, INF, (0.5, 0.0, 0.0), X_AXIS, 1.0)]
    kick = [Segment.rest(0.0, 1.0, ORIGIN), Segment.inertial(1.0, INF, ORIGIN, (0.5, 0.0, 0.0))]
    for broken in (gap, kick):
        w = Worldline(broken)  # building does not check: a query still has an answer
        assert not w.is_continuous()
        with pytest.raises(ValueError, match="do not join"):
            w.check_continuity()
    assert Worldline(gap).is_continuous(position_tol=0.6)
    assert Worldline(kick).is_continuous(velocity_tol=0.6)


def test_construction_errors():
    rest = Segment.rest(0.0, 1.0, ORIGIN)
    with pytest.raises(ValueError, match="at least one segment"):
        Worldline([])
    with pytest.raises(ValueError, match="order"):
        Worldline([Segment.rest(5.0, 6.0, ORIGIN), rest])
    with pytest.raises(ValueError, match="finite t0"):
        Worldline([Segment.rest(-INF, INF, ORIGIN)])
    with pytest.raises(ValueError):
        Worldline([rest], anchor_t=math.nan)
    with pytest.raises(TypeError, match="Segment"):
        Worldline([rest, "coast"])
    with pytest.raises(TypeError, match="Piece"):
        Worldline([plan.rest(plan.Start(), 1.0)])  # a plan's piece: pass its .segment
    with pytest.raises(TypeError):
        Worldline(rest)
    with pytest.raises(TypeError):
        Worldline()
    # a burn is not defined before it starts or after its mass is gone, so it cannot be an end segment
    burn = Segment.burn(1.0, 1.5, ORIGIN, X_AXIS, 0.0, 0.5, 0.8)
    for segments in ([burn], [rest, burn], [burn, Segment.rest(1.5, INF, ORIGIN)]):
        with pytest.raises(ValueError, match="burn"):
            Worldline(segments)
    end_x, end_v = burn.at(1.5)
    coast = Segment.inertial(1.5, INF, end_x, end_v)
    rocket = Worldline([Segment.rest(0.0, 1.0, ORIGIN), burn, coast])
    rocket.check_continuity()
    assert rocket.tau(1.5) == pytest.approx(1.0 + burn.dtau(1.5))
    # any iterable of segments will do
    assert len(Worldline(iter([rest]))) == 1 and len(Worldline((rest,))) == 1


# ------------------------------------------------------------------------------------------ vectorised


def expected(w: Worldline, ts, events):
    """The scalar answers, one call at a time."""
    xs = [w.at(t)[0] for t in ts]
    vs = [w.at(t)[1] for t in ts]
    taus = [w.tau(t) for t in ts]
    recvs = [w.recv(e[0], e[1:]) for e in events]
    sends = [w.send(e[0], e[1:]) for e in events]
    return xs, vs, taus, recvs, sends


def sample(w: Worldline, seed: int, n: int = 60):
    rng = random.Random(seed)
    ts = sample_times(w, rng, n)
    events = [(rng.uniform(-20.0, 120.0), *point(rng, 40.0)) for _ in range(n)]
    return ts, events


@pytest.mark.parametrize("make", [tour, proxima_ship])
def test_vectorised_methods_equal_the_scalar_ones_with_numpy(make):
    np = pytest.importorskip("numpy")
    w = make()
    ts, events = sample(w, 31)
    xs, vs, taus, recvs, sends = expected(w, ts, events)
    flat_events = [c for e in events for c in e]

    time_inputs = [
        ts,
        tuple(ts),
        array.array("d", ts),
        np.array(ts),
        np.array([ts, ts]).T.copy()[:, 1],  # a view that is not contiguous in memory
        np.array(ts)[::-1][::-1],
        memoryview(array.array("d", ts)),
        (t for t in ts),
    ]
    assert not time_inputs[4].flags["C_CONTIGUOUS"]
    for given in time_inputs:
        got_x, got_v = w.at_many(given)
        assert isinstance(got_x, np.ndarray) and got_x.dtype == np.float64
        assert got_x.shape == got_v.shape == (len(ts), 3)
        # element for element, bit for bit
        assert got_x.tolist() == [list(x) for x in xs]
        assert got_v.tolist() == [list(v) for v in vs]
    for given in time_inputs[:-1]:
        got = w.tau_many(given)
        assert isinstance(got, np.ndarray) and got.shape == (len(ts),)
        assert got.tolist() == taus

    event_inputs = [
        events,
        [list(e) for e in events],
        np.array(events),
        np.array(events).ravel(),
        flat_events,
        array.array("d", flat_events),
        np.asfortranarray(np.array(events)),
    ]
    for given in event_inputs:
        got_recv, got_send = w.recv_many(given), w.send_many(given)
        assert isinstance(got_recv, np.ndarray) and got_recv.shape == got_send.shape == (len(events),)
        assert got_recv.tolist() == recvs
        assert got_send.tolist() == sends


def test_vectorised_shapes_follow_the_input():
    np = pytest.importorskip("numpy")
    w = proxima_ship()
    grid = np.linspace(-1.0, 8.0, 24).reshape(2, 3, 4)
    xs, vs = w.at_many(grid)
    assert xs.shape == vs.shape == (2, 3, 4, 3)
    taus = w.tau_many(grid)
    assert taus.shape == (2, 3, 4)
    assert taus[1, 2, 3] == w.tau(float(grid[1, 2, 3]))
    assert tuple(xs[1, 0, 2]) == w.at(float(grid[1, 0, 2]))[0]
    events = np.zeros((5, 2, 4))
    events[..., 0] = np.linspace(0.0, 9.0, 10).reshape(5, 2)
    assert w.send_many(events).shape == (5, 2)
    assert w.send_many(events)[4, 1] == w.send(9.0, HOME)
    # a single time, an empty batch, whole numbers
    assert w.tau_many(2.0).shape == () and float(w.tau_many(2.0)) == w.tau(2.0)
    assert w.tau_many([]).shape == (0,) and w.at_many([])[0].shape == (0, 3)
    assert w.recv_many(np.zeros((0, 4))).shape == (0,)
    assert w.tau_many(range(4)).tolist() == [w.tau(float(t)) for t in range(4)]
    # the results are the caller's own: writable, and not shared between calls
    first, second = w.tau_many([1.0, 2.0]), w.tau_many([1.0, 2.0])
    first[0] = -1.0
    assert second[0] == w.tau(1.0)
    # an entry is infinite where the light never arrives
    runaway = Worldline([Segment.rest(-5.0, 0.0, ORIGIN), Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, 0.5)])
    assert runaway.recv_many([(0.25, -1.5, 0.0, 0.0), (0.75, -1.5, 0.0, 0.0)])[1] == INF

    for bad in (np.zeros((3, 3)), np.zeros(7), np.zeros(())):
        with pytest.raises(ValueError, match="four numbers"):
            w.recv_many(bad)
    with pytest.raises((TypeError, ValueError)):
        w.tau_many(["now"])


@pytest.mark.parametrize("make", [tour, proxima_ship])
def test_vectorised_methods_equal_the_scalar_ones_without_numpy(make, no_numpy):
    with pytest.raises(ImportError):
        import numpy  # noqa: F401  (the fixture hides it)

    w = make()
    ts, events = sample(w, 32)
    xs, vs, taus, recvs, sends = expected(w, ts, events)
    flat_events = [c for e in events for c in e]

    for given in (ts, tuple(ts), array.array("d", ts), (t for t in ts)):
        got_x, got_v = w.at_many(given)
        assert type(got_x) is list and type(got_v) is list
        # flat: three floats per time
        assert got_x == [c for x in xs for c in x]
        assert got_v == [c for v in vs for c in v]
    for given in (ts, array.array("d", ts), array.array("f", [1.0, 2.5])):
        got = w.tau_many(given)
        assert type(got) is list and all(type(c) is float for c in got)
        assert got == (taus if len(given) == len(ts) else [w.tau(1.0), w.tau(2.5)])

    for given in (events, [list(e) for e in events], flat_events, array.array("d", flat_events), iter(events)):
        assert w.recv_many(given) == recvs
    for given in (events, flat_events):
        assert w.send_many(given) == sends
    assert w.tau_many([]) == [] and w.at_many([]) == ([], []) and w.recv_many([]) == []
    assert w.tau_many(range(3)) == [w.tau(0.0), w.tau(1.0), w.tau(2.0)]

    with pytest.raises(ValueError, match="four numbers"):
        w.recv_many([(1.0, 2.0, 3.0)])
    with pytest.raises(ValueError, match="four numbers"):
        w.recv_many([1.0, 2.0, 3.0, 4.0, 5.0])
    with pytest.raises(TypeError):
        w.tau_many([(1.0, 2.0)])
    with pytest.raises(TypeError):
        w.tau_many(3.0)


def test_numpy_arrays_are_just_sequences_when_numpy_cannot_be_imported(monkeypatch):
    """An array made before NumPy went away is still a sequence of numbers; the answer is then a list."""
    import sys

    np = pytest.importorskip("numpy")
    ts = np.array([0.5, 1.0, 2.0])
    events = np.array([[3.0, 0.0, 0.0, 0.0], [5.0, 1.0, 0.0, 0.0]])
    monkeypatch.setitem(sys.modules, "numpy", None)
    w = proxima_ship()
    assert w.tau_many(ts) == [w.tau(0.5), w.tau(1.0), w.tau(2.0)]
    assert w.at_many(ts)[0][3:6] == list(w.at(1.0)[0])
    assert w.send_many(events) == [w.send(3.0, HOME), w.send(5.0, (1.0, 0.0, 0.0))]
