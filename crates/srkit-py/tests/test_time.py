"""Exact time: integer milliseconds to and from Python `int` without loss, and exact causal decisions."""

from __future__ import annotations

import math
import random

import pytest

import srkit
from srkit import Duration384, Duration512, Time384, Time512, causal

C = 299_792_458  # metres per second, exactly
PAIRS = [(Time384, Duration384, 384), (Time512, Duration512, 512)]


def interesting_ints(bits: int) -> list[int]:
    rng = random.Random(bits)
    values = [0, 1, 2, 999, 2**53, 2**63 - 1, 2**63, 2**64 - 1, 2**64, 2**64 + 1, 2**127, 2**128 - 1, 2**128]
    values += [10**k for k in (20, 50, 100, 115)]
    values += [2**bits - 1, 2**bits - 2, 2 ** (bits - 1), 2 ** (bits - 64) + 12345]
    values += [rng.getrandbits(rng.randrange(1, bits + 1)) for _ in range(300)]
    return values


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_python_ints_round_trip_without_loss(time, duration, bits):
    for n in interesting_ints(bits):
        t = time(n)
        assert int(t) == n and type(int(t)) is int
        assert str(t) == str(n) and repr(t) == f"{time.__name__}({n})"
        assert time.from_decimal(str(n)) == t and time(t) == t
        assert eval(repr(t), {time.__name__: time}) == t
        for signed in (n, -n):
            d = duration(signed)
            assert int(d) == signed and type(int(d)) is int
            assert str(d) == str(signed) and repr(d) == f"{duration.__name__}({signed})"
            assert duration.from_decimal(str(signed)) == d and duration(d) == d
            assert d.is_negative == (signed < 0) and bool(d) == (signed != 0)
    assert time.BITS == duration.BITS == bits
    assert int(time.MAX) == 2**bits - 1 and int(time.ZERO) == 0 and int(time()) == 0
    assert int(duration.ZERO) == 0 and int(duration()) == 0
    assert int(-duration(0)) == 0 and not duration(0).is_negative  # no negative zero


def test_the_examples_from_the_brief():
    assert int(Time512(10**150) + 1) == 10**150 + 1
    assert int(Time384(10**100) + 1) == 10**100 + 1
    # one light-second is exactly 1000 ticks, and one metre more is 1001
    assert causal.light_ticks((C, 0, 0)) == 1000
    assert causal.light_ticks((C + 1, 0, 0)) == 1001
    assert causal.C_METRES_PER_SECOND == C


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_arithmetic_is_exact_and_checked(time, duration, bits):
    rng = random.Random(bits + 1)
    top = 2**bits - 1
    for _ in range(300):
        a, b = rng.getrandbits(rng.randrange(1, bits + 1)), rng.getrandbits(rng.randrange(1, bits))
        ta, tb = time(a), time(b)
        # the difference of two times is a signed duration, exact across the whole width
        assert int(ta - tb) == a - b and isinstance(ta - tb, duration)
        assert int(tb - ta) == b - a
        assert ta - (ta - tb) == tb and tb + (ta - tb) == ta and (ta - tb) + tb == ta
        # a time plus or minus milliseconds, as an int or as a duration
        if a + b <= top:
            assert int(ta + b) == int(b + ta) == int(ta + duration(b)) == int(duration(b) + ta) == a + b
        else:
            with pytest.raises(OverflowError):
                ta + b
        if a >= b:
            assert int(ta - b) == int(ta - duration(b)) == int(ta + (-b)) == int(ta + duration(-b)) == a - b
        else:
            with pytest.raises(OverflowError):
                ta - b  # before the epoch
        # durations add and subtract like ints
        da, db = duration(a), duration(-b)
        assert int(da + db) == int(db + da) == a - b
        assert int(-da) == -a and int(abs(db)) == b and int(+db) == -b
        if a + b <= top:
            assert int(da - db) == a + b and int(db - da) == -(a + b)
        else:
            with pytest.raises(OverflowError):
                da - db
        if b + 5 <= top:
            assert int(5 - db) == 5 + b and int(db - 5) == -b - 5 and int(5 + db) == 5 - b
    with pytest.raises(OverflowError):
        time.MAX + 1
    with pytest.raises(OverflowError):
        time.ZERO - 1
    with pytest.raises(OverflowError):
        duration(top) + 1
    assert int(time.MAX - time.ZERO) == top and int(time.ZERO - time.MAX) == -top
    assert int(time.MAX - top) == 0


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_ordering_and_hashing(time, duration, bits):
    values = sorted(set(interesting_ints(bits)))
    times = [time(n) for n in values]
    assert sorted(reversed(times)) == times and max(times) == times[-1] and min(times) == times[0]
    for a, b in zip(times, times[1:]):
        assert a < b and a <= b and b > a and b >= a and a != b and not (a == b) and not (a > b)
    assert time(7) == time(7) and time(7) <= time(7) and hash(time(7)) == hash(time(7))
    assert len({time(n) for n in values} | {time(n) for n in values}) == len(values)
    assert {time(10**100): "late"}[time(10**100)] == "late"
    # a time is not a number: it does not equal the int it counts, and has no order against one
    assert time(7) != 7 and time(7) != "7" and time(7) != duration(7)
    with pytest.raises(TypeError):
        time(7) < 8

    signed = sorted({*values, *(-n for n in values)})
    durations = [duration(n) for n in signed]
    assert sorted(reversed(durations)) == durations
    for a, b in zip(durations, durations[1:]):
        assert a < b and b > a and a != b
    # a duration compares with an int of any size, and hashes as that int does
    assert duration(-7) == -7 and duration(7) > 6 and duration(7) < 2**600 and duration(7) > -(2**600)
    assert hash(duration(-7)) == hash(-7) and {duration(5): "x"}[5] == "x"
    assert duration(7) != 7.5 and duration(7) != "7"
    with pytest.raises(TypeError):
        duration(7) < 7.5


def test_the_two_widths_do_not_mix():
    assert Time384(5) != Time512(5) and Duration384(5) != Duration512(5)
    for mix in [
        lambda: Time384(5) - Time512(5),
        lambda: Time384(5) + Duration512(5),
        lambda: Duration384(5) + Duration512(5),
        lambda: Time384(5) < Time512(6),
        lambda: Time512(5).local_seconds_since(Time384(1)),
        lambda: causal.arrived(Time384(0), Time512(5), (0, 0, 0)),
        lambda: causal.settle_arrival(Time384(0), Time512(5), lambda tick: (0, 0, 0)),
        lambda: Time512(Time384(5)),
        lambda: Time384(5) + Time384(5),  # two times do not add
        lambda: 5 - Time384(5),
        lambda: Time384(5) * 2,
    ]:
        with pytest.raises(TypeError):
            mix()
    # the wider type holds what the narrower cannot
    assert int(Time512(2**384)) == 2**384
    with pytest.raises(OverflowError, match="384"):
        Time384(2**384)


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_what_is_not_an_exact_count_is_refused(time, duration, bits):
    with pytest.raises(OverflowError):
        time(2**bits)
    with pytest.raises(OverflowError):
        duration(-(2**bits))
    with pytest.raises(ValueError, match="negative"):
        time(-1)
    for not_an_int in (1.0, 1e30, "12", None, 1 + 0j):
        with pytest.raises(TypeError):
            time(not_an_int)
        with pytest.raises(TypeError):
            duration(not_an_int)
    with pytest.raises(TypeError):
        time(5) + 1.0
    for text in ("", "-1", "+1", " 1", "1 ", "1_000", "1e3", "0x10", "12.0"):
        with pytest.raises(ValueError):
            time.from_decimal(text)
    for text in ("", "--1", "+1", " 1", "1.5"):
        with pytest.raises(ValueError):
            duration.from_decimal(text)
    with pytest.raises(OverflowError):
        time.from_decimal(str(2**bits))
    assert int(duration.from_decimal("-0")) == 0 and int(time.from_decimal("000123")) == 123
    assert int(time(True)) == 1  # a bool is an int in Python


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_modular_phase_at_huge_epochs(time, duration, bits):
    rng = random.Random(bits + 2)
    year_ms = 31_557_600_000
    for _ in range(200):
        n = rng.getrandbits(rng.randrange(1, bits + 1))
        period = rng.choice([1, 2, 1000, year_ms, 86_400_000, rng.randrange(1, 2**64)])
        assert time(n) % period == n % period
    # an orbit's phase a googol of milliseconds on is exactly the phase of the remainder
    late = time(year_ms * 10**90 + 123_456_789)
    assert late % year_ms == 123_456_789
    assert late.phase_seconds(year_ms) == 123_456.789
    assert (late + year_ms) % year_ms == late % year_ms
    assert type(late % year_ms) is int
    with pytest.raises(ValueError):
        late % 0
    with pytest.raises(ValueError):
        late % -5
    with pytest.raises(OverflowError):
        late % 2**64
    with pytest.raises(TypeError):
        late % 1000.0
    with pytest.raises(ValueError):
        late.phase_seconds(0)


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_local_seconds_subtract_first(time, duration, bits):
    epoch = time(10**100)
    assert (epoch + 1500).local_seconds_since(epoch) == 1.5
    assert epoch.local_seconds_since(epoch + 250) == -0.25
    assert ((epoch + 1500) - epoch).local_seconds() == 1.5 and duration(-250).local_seconds() == -0.25
    # up to 2**53 ms every millisecond is a distinct float; beyond it the answer is refused, not rounded
    assert (epoch + 2**53).local_seconds_since(epoch) == 2**53 / 1000.0
    with pytest.raises(srkit.PrecisionLoss):
        (epoch + 2**53 + 1).local_seconds_since(epoch)
    with pytest.raises(srkit.PrecisionLoss):
        epoch.local_seconds_since(time(0))
    with pytest.raises(ArithmeticError):
        duration(10**30).local_seconds()
    assert issubclass(srkit.PrecisionLoss, ArithmeticError) and not issubclass(srkit.PrecisionLoss, ValueError)

    assert int((epoch + 5).duration_since(epoch)) == 5
    with pytest.raises(ValueError):
        epoch.duration_since(epoch + 5)  # the origin is later


# ----------------------------------------------------------------------------------------------- causal


def test_light_ticks_are_the_first_that_pass():
    assert causal.light_ticks((0, 0, 0)) == 0
    assert causal.light_ticks((1, 0, 0)) == 1  # any distance at all takes a tick
    assert causal.light_ticks((C, 0, 0)) == 1000 and causal.light_ticks((-C, 0, 0)) == 1000
    assert causal.light_ticks((C - 1, 0, 0)) == 1000 and causal.light_ticks((C + 1, 0, 0)) == 1001
    assert causal.light_ticks((3 * C, -4 * C, 0)) == 5000  # a 3-4-5 triangle of light-seconds
    assert type(causal.light_ticks((C, 0, 0))) is int
    rng = random.Random(51)
    for _ in range(500):
        scale = rng.randrange(0, 60)
        dx = tuple(rng.randrange(-(2**63), 2**63) << scale for _ in range(3))
        ticks = causal.light_ticks(dx)
        # the smallest n with (c n)^2 >= 10^6 |dx|^2, by Python's own integers
        norm2 = 10**6 * sum(c * c for c in dx)
        assert (C * ticks) ** 2 >= norm2 and (ticks == 0 or (C * (ticks - 1)) ** 2 < norm2)
        # along one axis it is the plain ceiling
        assert causal.light_ticks((dx[0], 0, 0)) == -((-abs(dx[0]) * 1000) // C)
    # the largest displacement there is: nothing overflows
    far = (2**127 - 1, -(2**127), 2**127 - 1)
    ticks = causal.light_ticks(far)
    assert (C * ticks) ** 2 >= 10**6 * sum(c * c for c in far) > (C * (ticks - 1)) ** 2


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_arrived_flips_exactly_at_first_arrival(time, duration, bits):
    rng = random.Random(bits + 3)
    epochs = [0, 12_345, 10**12, 10**100, 7 * 10**110, 2 ** (bits - 1)]
    for _ in range(100):
        dx = tuple(rng.randrange(-(10**27), 10**27) >> rng.randrange(0, 80) for _ in range(3))
        ticks = causal.light_ticks(dx)
        for n in epochs:
            emission = time(n)
            arrival = causal.first_arrival(emission, dx)
            assert type(arrival) is time and arrival - emission == ticks  # the same at any epoch
            assert causal.arrived(emission, arrival, dx)
            assert causal.arrived(emission, arrival + 1, dx) and causal.arrived(emission, arrival + 10**9, dx)
            if ticks:
                assert not causal.arrived(emission, arrival - 1, dx)
            # before the emission nothing has arrived, even at no distance
            if n:
                assert not causal.arrived(emission, emission - 1, (0, 0, 0))
            assert causal.arrived(emission, emission, (0, 0, 0))
            # the stationary-distance forms on the time types agree along an axis
            distance = abs(dx[0])
            assert emission.stationary_arrival(distance) == causal.first_arrival(emission, (dx[0], 0, 0))
            assert int(duration.light_travel_time(distance)) == causal.light_ticks((dx[0], 0, 0))
            late = emission.stationary_arrival(distance)
            assert late.stationary_signal_arrived(emission, distance)
            assert not (distance and (late - 1).stationary_signal_arrived(emission, distance))
    # the far future, across the observable universe (8.8e26 m)
    emission, far = Time512(7 * 10**140), (880 * 10**24, 0, 0)
    arrival = causal.first_arrival(emission, far)
    assert causal.arrived(emission, arrival, far) and not causal.arrived(emission, arrival - 1, far)
    with pytest.raises(OverflowError):
        causal.first_arrival(time.MAX, (C, 0, 0))


def test_classification_is_exact_on_the_cone():
    one_second = (C, 0, 0)
    # c * 1000 ms = 1 light-second exactly: null, and one metre either way tips it
    assert causal.classify(1000, one_second) == "null"
    assert causal.classify(1000, (C - 1, 0, 0)) == "timelike"
    assert causal.classify(1000, (C, 1, 0)) == "spacelike"
    assert causal.classify(999, one_second) == "spacelike" and causal.classify(1001, one_second) == "timelike"
    # the sign of the interval does not matter; durations of either width are accepted
    assert causal.classify(-1000, (0, 0, -C)) == "null"
    assert causal.classify(Duration384(-1000), one_second) == "null"
    assert causal.classify(Duration512(1000), one_second) == "null"
    assert causal.classify(Time384(10**100 + 1000) - Time384(10**100), one_second) == "null"
    assert causal.classify(0, (0, 0, 0)) == "null" and causal.classify(0, (1, 0, 0)) == "spacelike"
    # a wait far longer than any distance is timelike without overflowing
    assert causal.classify(2**512 - 1, (2**127 - 1, -(2**127), 2**127 - 1)) == "timelike"
    # a 3-4-5 triangle scaled up towards the edge of the grid: still exact
    k = 10**20
    assert causal.classify(5000 * k, (3 * C * k, 4 * C * k, 0)) == "null"
    assert causal.classify(5000 * k, (3 * C * k, 4 * C * k, 1)) == "spacelike"
    assert causal.classify(5000 * k + 1, (3 * C * k, 4 * C * k, 1)) == "timelike"


def test_numpy_integers_are_exact_and_numpy_floats_are_not():
    np = pytest.importorskip("numpy")
    assert causal.light_ticks(np.array([C, 0, 0])) == 1000
    assert causal.light_ticks((np.int64(C), np.int32(1), np.uint8(0))) == 1001
    assert causal.classify(np.int64(1000), np.array([C, 0, 0], dtype=np.int64)) == "null"
    for not_exact in [np.array([1.0, 0.0, 0.0]), (np.float64(1.0), 0, 0)]:
        with pytest.raises(TypeError):
            causal.light_ticks(not_exact)
    with pytest.raises(TypeError):
        Time384(np.float64(5.0))


def test_displacements_are_whole_metres():
    assert causal.light_ticks([C, 0, 0]) == causal.light_ticks((C, 0, 0)) == 1000
    assert causal.light_ticks((C, True, False)) == 1001  # (C, 1, 0): one metre off the axis
    for not_exact in [(1.0, 0, 0), (C, 0.5, 0), ("1", 0, 0), (None, 0, 0)]:
        with pytest.raises(TypeError):
            causal.light_ticks(not_exact)
    with pytest.raises(TypeError):
        causal.light_ticks(5)
    with pytest.raises(TypeError):
        causal.light_ticks("abc")
    with pytest.raises(ValueError, match="three"):
        causal.light_ticks((1, 2))
    # each component must fit a signed 128-bit integer of metres
    assert causal.light_ticks((2**127 - 1, 0, 0)) > 0 and causal.light_ticks((-(2**127), 0, 0)) > 0
    for too_far in [(2**127, 0, 0), (0, -(2**127) - 1, 0), (0, 0, 10**60)]:
        with pytest.raises(OverflowError):
            causal.light_ticks(too_far)
        with pytest.raises(OverflowError):
            causal.classify(1, too_far)
    with pytest.raises(OverflowError):
        causal.classify(2**512, (0, 0, 0))
    with pytest.raises(TypeError):
        causal.classify(1.0, (0, 0, 0))
    with pytest.raises(TypeError, match="Time"):
        causal.first_arrival(10**100, (C, 0, 0))  # an int is not a time: say which width

    assert causal.grid_from_metres(1.4, -2.6, 0.5) == (1, -3, 1)  # to nearest, ties away from zero
    assert causal.grid_from_metres(-0.5, 1e30, 0.0) == (-1, int(1e30), 0)
    assert all(type(c) is int for c in causal.grid_from_metres(1.0, 2.0, 3.0))
    for bad in [(math.nan, 0.0, 0.0), (math.inf, 0.0, 0.0), (1e39, 0.0, 0.0)]:
        with pytest.raises(OverflowError):
            causal.grid_from_metres(*bad)


@pytest.mark.parametrize("time, duration, bits", PAIRS)
def test_settling_an_arrival_on_a_moving_receiver(time, duration, bits):
    emission = time(10**60)
    c_m_per_ms = C / 1000.0
    v = 0.6 * c_m_per_ms  # the receiver recedes at 0.6 c, metres per millisecond
    x0 = 4.0e9  # metres away at emission

    def displacement(tick):
        dt_ms = tick.local_seconds_since(emission) * 1e3
        return causal.grid_from_metres(x0 + v * dt_ms, 0.0, 0.0)

    exact_ms = x0 / (c_m_per_ms - v)
    answers = set()
    for off in range(-6, 7):
        estimate = emission + (math.ceil(exact_ms) + off)
        settled = causal.settle_arrival(emission, estimate, displacement)
        assert type(settled) is time
        assert causal.arrived(emission, settled, displacement(settled))
        assert not causal.arrived(emission, settled - 1, displacement(settled - 1))
        answers.add(int(settled))
    assert len(answers) == 1  # the same first tick from any estimate within the guard
    assert abs(answers.pop() - int(emission) - exact_ms) <= 1.0

    # an estimate further off than the guard is refused, not searched for ever
    with pytest.raises(srkit.PrecisionLoss):
        causal.settle_arrival(emission, emission + (math.ceil(exact_ms) + 500), displacement)
    assert causal.settle_arrival(emission, emission + (math.ceil(exact_ms) + 500), displacement, guard=1000)

    # the callback's own errors come through as they are
    def broken(tick):
        raise RuntimeError("no ephemeris")

    with pytest.raises(RuntimeError, match="no ephemeris"):
        causal.settle_arrival(emission, emission + 5, broken)
    with pytest.raises(TypeError):
        causal.settle_arrival(emission, emission + 5, lambda tick: (1.5, 0, 0))
    with pytest.raises(TypeError):
        causal.settle_arrival(emission, emission + 5, "not callable")
    # an estimate before the emission starts from the emission
    assert causal.settle_arrival(emission, emission - 5, lambda tick: (0, 0, 0)) == emission
