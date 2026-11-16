"""`srkit.causal.Clock`: a width, a tick and a grid spacing, and the exact checks at that resolution."""

from __future__ import annotations

import math
import random

import pytest

import srkit
from srkit import (
    Duration384,
    Duration512,
    Duration768,
    Duration1024,
    Time384,
    Time512,
    Time768,
    Time1024,
    causal,
)
from srkit.causal import Clock

C = 299_792_458  # metres per second, exactly
WIDTHS = {
    384: (Time384, Duration384),
    512: (Time512, Duration512),
    768: (Time768, Duration768),
    1024: (Time1024, Duration1024),
}
# (ticks per second, grid steps per metre): the three named ones, the finest allowed, and some odd ones
RESOLUTIONS = [(1000, 1), (10**6, 1), (10**9, 1000), (10**12, 10**6), (1, 1), (60, 100), (1024, 3)]
BITS = pytest.mark.parametrize("bits", sorted(WIDTHS))


def first_tick(ticks_per_second: int, grid_per_metre: int, dx) -> int:
    """The smallest n with c g n >= f |dx|, by Python's own integers."""
    norm2 = ticks_per_second**2 * sum(c * c for c in dx)
    n = math.isqrt(norm2) // (C * grid_per_metre)
    while (C * grid_per_metre * n) ** 2 < norm2:
        n += 1
    return n


# --------------------------------------------------------------------------------------- construction


def test_the_default_clock_is_384_bits_of_milliseconds_on_a_metre_grid():
    clock = Clock()
    assert (clock.bits, clock.ticks_per_second, clock.grid_per_metre) == (384, 1000, 1)
    assert clock == Clock(384) == Clock(384, 1000, 1) == Clock(bits=384, ticks_per_second=1000, grid_per_metre=1)
    assert clock == Clock(384, *causal.MILLISECOND_METRE)
    assert causal.MILLISECOND_METRE == (1000, 1) and type(causal.MILLISECOND_METRE) is tuple
    assert causal.MICROSECOND_METRE == (10**6, 1)
    assert causal.NANOSECOND_MILLIMETRE == (10**9, 1000)
    assert all(type(n) is int for pair in (causal.MILLISECOND_METRE, causal.NANOSECOND_MILLIMETRE) for n in pair)


@BITS
def test_properties(bits):
    time, duration = WIDTHS[bits]
    for ticks, grid in RESOLUTIONS:
        clock = Clock(bits, ticks, grid)
        assert (clock.bits, clock.ticks_per_second, clock.grid_per_metre) == (bits, ticks, grid)
        assert all(type(n) is int for n in (clock.bits, clock.ticks_per_second, clock.grid_per_metre))
        assert clock.tick_seconds == 1.0 / ticks and clock.grid_metres == 1.0 / grid
        assert clock.light_steps_per_tick == pytest.approx(C * grid / ticks, rel=1e-15)
        # the longest time it holds: 2**bits ticks, in Julian years
        years = bits * math.log10(2) - math.log10(ticks) - math.log10(31_557_600)
        assert clock.reach_log10_years == pytest.approx(years, abs=1e-9)
        assert clock.Time is time and clock.Duration is duration
        assert repr(clock) == f"Clock(bits={bits}, ticks_per_second={ticks}, grid_per_metre={grid})"
        assert eval(repr(clock), {"Clock": Clock}) == clock
    with pytest.raises(AttributeError):
        Clock(bits).bits = 512  # a clock does not change


def test_reach():
    assert Clock(768, 10**9, 1000).reach_log10_years == pytest.approx(214.7, abs=0.05)
    assert Clock(384).reach_log10_years == pytest.approx(105.1, abs=0.05)
    # 384 bits of nanoseconds fall just short of a googol of years; 512 bits do not
    assert Clock(384, 10**9, 1000).reach_log10_years < 100 < Clock(512, 10**9, 1000).reach_log10_years
    assert Clock(1024, *causal.NANOSECOND_MILLIMETRE).reach_log10_years == pytest.approx(291.76, abs=0.01)
    # and the reach is the truth: MAX is that many years of ticks
    for bits, (time, _) in WIDTHS.items():
        clock = Clock(bits, 10**6, 1)
        years = int(time.MAX) // (10**6 * 31_557_600)
        assert math.log10(years) == pytest.approx(clock.reach_log10_years, abs=1e-9)


def test_clocks_compare_and_hash_by_their_three_numbers():
    a, b = Clock(768, 10**9, 1000), Clock(768, *causal.NANOSECOND_MILLIMETRE)
    assert a == b and not (a != b) and hash(a) == hash(b) and a is not b
    assert len({a, b, Clock(768, 10**9, 1000)}) == 1 and {a: "x"}[b] == "x"
    others = [Clock(384, 10**9, 1000), Clock(768, 10**6, 1000), Clock(768, 10**9, 1), Clock(1024, 10**9, 1000)]
    assert all(a != other and not (a == other) for other in others)
    assert len({a, *others}) == 5
    assert a != (768, 10**9, 1000) and a != "clock" and a != 768 and a is not None and a != None  # noqa: E711
    with pytest.raises(TypeError):
        a < b


def test_bad_bits_and_resolutions_are_refused():
    for bits in (0, 1, 64, 128, 256, 383, 385, 640, 2048, -384, 2**70):
        with pytest.raises(ValueError, match="384, 512, 768 or 1024"):
            Clock(bits)
    for ticks, grid in [(0, 1), (1000, 0), (-1000, 1), (1000, -1), (10**12 + 1, 1), (1000, 10**6 + 1), (2**64, 1), (1, 2**70)]:
        with pytest.raises(ValueError, match="ticks_per_second"):
            Clock(384, ticks, grid)
    Clock(1024, 10**12, 10**6)  # the finest allowed: a picosecond and a micrometre
    Clock(384, 1, 1)  # and the coarsest
    # the three numbers are exact counts: floats, text and None are not
    for bad in [(384.0,), ("384",), (None,), (384, 1000.0), (384, 1e9, 1000), (384, 1000, 1.0), (384, None), (384, 1000, None)]:
        with pytest.raises(TypeError, match="must be an int"):
            Clock(*bad)
    with pytest.raises(TypeError):
        Clock(384, 1000, 1, 1)
    assert Clock(True + 383).bits == 384


def test_numpy_integers_are_counts_too():
    np = pytest.importorskip("numpy")
    clock = Clock(np.int64(768), np.int64(10**9), np.uint16(1000))
    assert clock == Clock(768, 10**9, 1000) and type(clock.bits) is int
    assert clock.light_ticks(np.array([C * 1000, 0, 0])) == 10**9
    assert int(clock.time(np.int64(5))) == 5 and int(clock.seconds(np.int32(-2))) == -2 * 10**9
    with pytest.raises(TypeError):
        Clock(np.float64(768.0))
    with pytest.raises(TypeError):
        clock.light_ticks(np.array([1.0, 0.0, 0.0]))


# ------------------------------------------------------------------------------- times and durations


@BITS
def test_times_and_durations_of_every_width(bits):
    time, duration = WIDTHS[bits]
    top = 2**bits - 1
    assert time.BITS == duration.BITS == bits and int(time.MAX) == top
    rng = random.Random(bits)
    values = [0, 1, 2**64 - 1, 2**64, 2 ** (bits - 1), top - 1, top]
    values += [rng.getrandbits(rng.randrange(1, bits + 1)) for _ in range(100)]
    for ticks, grid in [(1000, 1), (10**9, 1000), (10**12, 10**6)]:
        clock = Clock(bits, ticks, grid)
        for n in values:
            # a time is a count of ticks whatever their length: the clock builds the class of its width
            t = clock.time(n)
            assert type(t) is time and int(t) == n and t == time(n) and type(int(t)) is int
            for signed in (n, -n):
                d = clock.duration(signed)
                assert type(d) is duration and int(d) == signed and d == duration(signed)
        assert clock.time(top) == time.MAX
        with pytest.raises(OverflowError):
            clock.time(top + 1)
        with pytest.raises(OverflowError):
            clock.time(top) + 1
        with pytest.raises(OverflowError):
            clock.duration(-(top + 1))
        with pytest.raises(ValueError, match="negative"):
            clock.time(-1)
        for not_a_count in (1.0, "5", None, time(5), duration(5)):
            with pytest.raises(TypeError):
                clock.time(not_a_count)
            with pytest.raises(TypeError):
                clock.duration(not_a_count)


@BITS
def test_whole_seconds(bits):
    _, duration = WIDTHS[bits]
    for ticks, grid in RESOLUTIONS:
        clock = Clock(bits, ticks, grid)
        assert type(clock.seconds(1)) is duration and clock.seconds(1) == ticks
        assert clock.seconds(0) == 0 and not clock.seconds(0).is_negative
        assert clock.seconds(-90) == -90 * ticks and clock.seconds(31_557_600) == 31_557_600 * ticks
        assert clock.local_seconds(clock.seconds(3600)) == 3600.0
        # far more seconds than 64 bits count: a googol of years, where the width holds it
        googol_years = 10**100 * 31_557_600
        if googol_years * ticks < 2**bits:
            assert int(clock.seconds(googol_years)) == googol_years * ticks
        else:
            with pytest.raises(OverflowError):
                clock.seconds(googol_years)
        most = (2**bits - 1) // ticks
        assert int(clock.seconds(most)) == most * ticks and int(clock.seconds(-most)) == -most * ticks
        if ticks > 1:
            with pytest.raises(OverflowError):
                clock.seconds(most + 1)
        with pytest.raises(TypeError):
            clock.seconds(1.5)


@BITS
def test_local_seconds_subtract_first(bits):
    time, duration = WIDTHS[bits]
    for ticks, grid in RESOLUTIONS:
        clock = Clock(bits, ticks, grid)
        epoch = clock.time(10**100)
        assert clock.local_seconds(clock.seconds(3)) == 3.0 and clock.local_seconds(clock.seconds(-3)) == -3.0
        assert clock.local_seconds(7) == 7 / ticks and clock.local_seconds(duration(-7)) == -7 / ticks
        assert clock.local_seconds_between(epoch + 3 * ticks, epoch) == 3.0
        assert clock.local_seconds_between(epoch, epoch + 3 * ticks) == -3.0
        assert clock.local_seconds_between(epoch + 1, epoch) == 1 / ticks == clock.tick_seconds
        # up to 2**53 ticks every tick is a distinct float; beyond it the answer is refused, not rounded
        assert clock.local_seconds(2**53) == 2**53 / ticks
        assert clock.local_seconds(-(2**53)) == -(2**53) / ticks
        assert clock.local_seconds_between(epoch + 2**53, epoch) == 2**53 / ticks
        for too_long in (2**53 + 1, -(2**53) - 1, 10**30, duration(2**53 + 1)):
            with pytest.raises(srkit.PrecisionLoss):
                clock.local_seconds(too_long)
        with pytest.raises(srkit.PrecisionLoss):
            clock.local_seconds_between(epoch + (2**53 + 1), epoch)
        with pytest.raises(srkit.PrecisionLoss):
            clock.local_seconds_between(epoch, time(0))
        with pytest.raises(TypeError):
            clock.local_seconds(1.5)
        with pytest.raises(TypeError):
            clock.local_seconds_between(10**100 + 5, epoch)  # an int is not a time
    # at the default resolution the methods on the classes say the same
    clock, epoch = Clock(bits), time(10**100)
    assert clock.local_seconds_between(epoch + 1500, epoch) == (epoch + 1500).local_seconds_since(epoch) == 1.5
    assert clock.local_seconds(duration(-250)) == duration(-250).local_seconds() == -0.25


@BITS
def test_phase_at_huge_epochs(bits):
    time, _ = WIDTHS[bits]
    rng = random.Random(bits + 7)
    for ticks, grid in RESOLUTIONS:
        clock = Clock(bits, ticks, grid)
        day = 86_400 * ticks  # a day of ticks: the period is a 64-bit count, the epoch is not
        late = clock.time(day * 10**80 + 123_456_789)
        assert clock.phase_seconds(late, day) == (123_456_789 % day) / ticks
        assert clock.phase_seconds(late + day, day) == clock.phase_seconds(late, day)
        for _ in range(50):
            n, period = rng.getrandbits(rng.randrange(1, bits + 1)), rng.randrange(1, 2**53)
            assert clock.phase_seconds(time(n), period) == (n % period) / ticks
        for bad_period in (0, -5):
            with pytest.raises(ValueError):
                clock.phase_seconds(late, bad_period)
        with pytest.raises(OverflowError):
            clock.phase_seconds(late, 2**64)
        with pytest.raises(TypeError):
            clock.phase_seconds(late, 1000.0)
        # a phase beyond 2**53 ticks cannot be a float to the tick
        with pytest.raises(srkit.PrecisionLoss):
            clock.phase_seconds(time(2**60), 2**61)
    assert Clock(bits).phase_seconds(time(10**90 + 1500), 1000) == time(10**90 + 1500).phase_seconds(1000) == 0.5


# ---------------------------------------------------------------------------------- exact decisions


@BITS
@pytest.mark.parametrize("ticks, grid", RESOLUTIONS)
def test_one_light_second_is_exactly_one_second_of_ticks(bits, ticks, grid):
    time, duration = WIDTHS[bits]
    clock = Clock(bits, ticks, grid)
    away = (C * grid, 0, 0)  # one light-second, in grid steps
    assert clock.light_ticks(away) == ticks and type(clock.light_ticks(away)) is int
    assert clock.seconds(1) == clock.light_ticks(away)
    assert clock.grid(float(C), 0.0, 0.0) == away
    for sent in (clock.time(0), clock.time(12_345), clock.time(3 * 10**39), clock.time(2 ** (bits - 1))):
        heard = clock.first_arrival(sent, away)
        assert type(heard) is time and heard - sent == ticks and type(heard - sent) is duration
        # the tick before, the signal has not arrived; from the exact tick on, it has
        assert not clock.arrived(sent, heard - 1, away)
        assert clock.arrived(sent, heard, away) and clock.arrived(sent, heard + 1, away)
        assert clock.arrived(sent, sent, (0, 0, 0)) and not clock.arrived(heard, sent, (0, 0, 0))
        # exactly on the light cone, and a tick either side
        assert clock.classify(heard - sent, away) == "null"
        assert clock.classify(heard - 1 - sent, away) == "spacelike"
        assert clock.classify(heard + 1 - sent, away) == "timelike"
    assert clock.classify(ticks, away) == clock.classify(-ticks, away) == "null"
    assert clock.classify(ticks - 1, away) == "spacelike" and clock.classify(ticks + 1, away) == "timelike"
    # one grid step more or less tips it
    assert clock.classify(ticks, (C * grid, 1, 0)) == "spacelike"
    assert clock.classify(ticks, (C * grid - 1, 0, 0)) == "timelike"
    assert clock.light_ticks((C * grid + 1, 0, 0)) == ticks + 1
    # a 3-4-5 triangle of light-seconds, in any direction
    assert clock.light_ticks((3 * C * grid, 0, -4 * C * grid)) == 5 * ticks
    assert clock.classify(5 * ticks, (0, -3 * C * grid, 4 * C * grid)) == "null"
    with pytest.raises(OverflowError):
        clock.first_arrival(time.MAX, away)


@pytest.mark.parametrize("ticks, grid", RESOLUTIONS)
def test_the_widths_agree(ticks, grid):
    rng = random.Random(ticks + grid)
    clocks = [Clock(bits, ticks, grid) for bits in sorted(WIDTHS)]
    for _ in range(200):
        dx = tuple(rng.randrange(-(2**126), 2**126) >> rng.randrange(0, 120) for _ in range(3))
        expected = first_tick(ticks, grid, dx)
        epoch = rng.getrandbits(rng.randrange(1, 380))
        dt = max(0, expected + rng.choice([-2, -1, 0, 1, 2]))
        relation = "null" if (C * grid * dt) ** 2 == ticks**2 * sum(c * c for c in dx) else None
        relation = relation or ("timelike" if dt >= expected else "spacelike")
        for clock in clocks:
            assert clock.light_ticks(dx) == expected
            sent = clock.time(epoch)
            heard = clock.first_arrival(sent, dx)
            assert int(heard) == epoch + expected
            assert clock.arrived(sent, sent + dt, dx) == (dt >= expected)
            assert clock.classify(dt, dx) == clock.classify(clock.duration(-dt), dx) == relation
    # the largest displacement there is: nothing overflows, at any resolution
    far = (2**127 - 1, -(2**127), 2**127 - 1)
    expected = first_tick(ticks, grid, far)
    for clock in clocks:
        assert clock.light_ticks(far) == expected
        assert clock.classify(expected, far) in ("timelike", "null") and clock.classify(expected - 1, far) == "spacelike"
        assert clock.classify(2**clock.bits - 1, far) == "timelike"
        with pytest.raises(OverflowError):
            clock.classify(2**clock.bits, far)


@BITS
def test_the_module_functions_are_the_default_clock(bits):
    time, duration = WIDTHS[bits]
    clock = Clock(bits)
    rng = random.Random(bits + 11)
    assert Clock(bits, 1000, 1) == clock == Clock(bits, *causal.MILLISECOND_METRE)
    for _ in range(200):
        dx = tuple(rng.randrange(-(10**27), 10**27) >> rng.randrange(0, 80) for _ in range(3))
        sent = time(rng.getrandbits(rng.randrange(1, bits)))
        ticks = causal.light_ticks(dx)
        assert clock.light_ticks(dx) == ticks and type(clock.light_ticks(dx)) is type(ticks) is int
        heard = causal.first_arrival(sent, dx)
        assert clock.first_arrival(sent, dx) == heard and type(heard) is time
        for now in (heard - 1 if ticks else heard, heard, heard + 1):
            assert clock.arrived(sent, now, dx) == causal.arrived(sent, now, dx)
            assert clock.classify(now - sent, dx) == causal.classify(now - sent, dx)
            assert clock.classify(int(now - sent), dx) == causal.classify(int(now - sent), dx)
        metres = tuple(float(c) for c in dx)
        assert clock.grid(*metres) == causal.grid_from_metres(*metres)

    def displacement(tick):
        return (4_000_000_000 + 180_000 * int(tick - sent), 0, 0)  # receding at 0.6 c, near enough

    for off in (-3, 0, 4):
        estimate = sent + (33_357 + off)
        settled = causal.settle_arrival(sent, estimate, displacement, guard=64)
        assert clock.settle_arrival(sent, estimate, displacement, guard=64) == settled


def test_grid_rounds_metres_to_grid_steps():
    assert Clock().grid(1.4, -2.6, 0.5) == (1, -3, 1) == causal.grid_from_metres(1.4, -2.6, 0.5)
    millimetres = Clock(768, *causal.NANOSECOND_MILLIMETRE)
    assert millimetres.grid(1.0, -0.0014, 0.0005) == (1000, -1, 1)  # to nearest, ties away from zero
    assert millimetres.grid(-0.0005, 0.00049, 2.5) == (-1, 0, 2500)
    assert Clock(384, 60, 100).grid(1.234, 0.0, -7.0) == (123, 0, -700)
    assert Clock(384, 1024, 3).grid(1.0, 0.5, 0.49) == (3, 2, 1)
    assert all(type(c) is int for c in millimetres.grid(1.0, 2.0, 3.0))
    # a light-year in micrometres is still an int, though no longer exact to the step
    far = Clock(1024, 10**12, 10**6).grid(9.4607304725808e15, 0.0, 0.0)
    assert far[0] == pytest.approx(9_460_730_472_580_800 * 10**6, rel=1e-15) and type(far[0]) is int
    for bad in [(math.nan, 0.0, 0.0), (0.0, math.inf, 0.0), (0.0, 0.0, 1e36)]:
        with pytest.raises(OverflowError):
            millimetres.grid(*bad)
    assert Clock().grid(1e36, 0.0, 0.0) == (int(1e36), 0, 0)  # the same metres fit the metre grid
    with pytest.raises(TypeError):
        millimetres.grid("1", 0.0, 0.0)


def test_displacements_are_whole_grid_steps():
    clock = Clock(512, 10**6, 1)
    assert clock.light_ticks([C, 0, 0]) == clock.light_ticks((C, 0, 0)) == 10**6
    for not_exact in [(1.0, 0, 0), (C, 0.5, 0), ("1", 0, 0), (None, 0, 0), 5, "abc"]:
        with pytest.raises(TypeError):
            clock.light_ticks(not_exact)
        with pytest.raises(TypeError):
            clock.classify(1, not_exact)
        with pytest.raises(TypeError):
            clock.first_arrival(clock.time(0), not_exact)
    with pytest.raises(ValueError, match="three"):
        clock.light_ticks((1, 2))
    for too_far in [(2**127, 0, 0), (0, -(2**127) - 1, 0)]:
        with pytest.raises(OverflowError):
            clock.light_ticks(too_far)
    with pytest.raises(TypeError):
        clock.classify(1.0, (0, 0, 0))


def test_a_clock_takes_only_times_of_its_own_width():
    for bits, (time, duration) in WIDTHS.items():
        clock = Clock(bits, 10**6, 1)
        for other_bits, (other_time, other_duration) in WIDTHS.items():
            if other_bits == bits:
                continue
            t, d = other_time(5), other_duration(5)
            for mix in [
                lambda: clock.first_arrival(t, (0, 0, 0)),
                lambda: clock.arrived(t, t, (0, 0, 0)),
                lambda: clock.arrived(time(0), t, (0, 0, 0)),
                lambda: clock.arrived(t, time(9), (0, 0, 0)),
                lambda: clock.classify(d, (0, 0, 0)),
                lambda: clock.settle_arrival(t, t, lambda tick: (0, 0, 0)),
                lambda: clock.settle_arrival(time(0), t, lambda tick: (0, 0, 0)),
                lambda: clock.local_seconds(d),
                lambda: clock.local_seconds_between(t, t),
                lambda: clock.local_seconds_between(time(9), t),
                lambda: clock.phase_seconds(t, 1000),
                lambda: clock.time(t),
                # and the classes and the module-level functions, as before
                lambda: time(5) - t,
                lambda: time(5) + d,
                lambda: duration(5) + d,
                lambda: time(5) < t,
                lambda: time(t),
                lambda: duration(d),
                lambda: causal.arrived(time(0), t, (0, 0, 0)),
                lambda: causal.settle_arrival(time(0), t, lambda tick: (0, 0, 0)),
            ]:
                with pytest.raises(TypeError):
                    mix()
            assert time(5) != t and duration(5) != d
        with pytest.raises(TypeError, match=time.__name__):
            clock.first_arrival(10**100, (C, 0, 0))  # an int is not a time: the message says which class
        with pytest.raises(TypeError, match=time.__name__):
            clock.arrived(time(0), Time384(5) if bits != 384 else Time512(5), (0, 0, 0))
    # each wider class holds what the narrower cannot
    assert int(Time768(2**512)) == 2**512 and int(Time1024(2**768)) == 2**768
    with pytest.raises(OverflowError, match="768"):
        Time768(2**768)
    with pytest.raises(OverflowError, match="1024"):
        Time1024(2**1024)
    with pytest.raises(TypeError, match="Time1024"):
        causal.first_arrival(10**100, (C, 0, 0))  # the module-level message lists every width


@BITS
@pytest.mark.parametrize("ticks, grid", [(1000, 1), (10**6, 1), (10**9, 1000), (10**12, 10**6), (60, 100)])
def test_settling_an_arrival_on_a_moving_receiver(bits, ticks, grid):
    time, _ = WIDTHS[bits]
    clock = Clock(bits, ticks, grid)
    emission = clock.time(10**60)
    x0 = 4.0e6  # metres away at the emission
    v = 0.6 * C  # receding, metres per second

    def displacement(tick):
        return clock.grid(x0 + v * clock.local_seconds_between(tick, emission), 0.0, 0.0)

    def arrived_by(n: int) -> bool:
        """The exact check by Python's own integers, n ticks after the emission."""
        dx = displacement(emission + n)[0]
        return C * grid * n >= ticks * dx

    exact = x0 / (C - v) * ticks  # the float estimate, in ticks
    first = math.floor(exact) - 4
    assert not arrived_by(first)
    while not arrived_by(first):
        first += 1
    assert abs(first - exact) <= 2.0

    for off in (-5, -2, -1, 0, 1, 3, 6):
        settled = clock.settle_arrival(emission, emission + (math.ceil(exact) + off), displacement)
        assert type(settled) is time and settled - emission == first  # the exact tick, from any estimate nearby
        assert clock.arrived(emission, settled, displacement(settled))
        assert not clock.arrived(emission, settled - 1, displacement(settled - 1))

    # an estimate further off than the guard is refused, not searched for ever
    with pytest.raises(srkit.PrecisionLoss):
        clock.settle_arrival(emission, emission + (first + 500), displacement)
    with pytest.raises(srkit.PrecisionLoss):
        clock.settle_arrival(emission, emission + (first + 3), displacement, guard=2)
    assert clock.settle_arrival(emission, emission + (first + 500), displacement, guard=1000) - emission == first
    assert clock.settle_arrival(emission, emission + (first + 3), displacement, guard=4) - emission == first
    assert clock.settle_arrival(emission, estimate=emission + first, displacement_at=displacement, guard=1) - emission == first

    # the callback is given times of the clock's width, and its own errors come through as they are
    seen = []
    clock.settle_arrival(emission, emission + first, lambda tick: seen.append(tick) or displacement(tick))
    assert seen and all(type(tick) is time for tick in seen)

    def broken(tick):
        raise RuntimeError("no ephemeris")

    with pytest.raises(RuntimeError, match="no ephemeris"):
        clock.settle_arrival(emission, emission + 5, broken)
    with pytest.raises(TypeError):
        clock.settle_arrival(emission, emission + 5, lambda tick: (1.5, 0, 0))
    with pytest.raises(TypeError):
        clock.settle_arrival(emission, emission + 5, "not callable")
    # an estimate before the emission starts from the emission
    assert clock.settle_arrival(emission, emission - 5, lambda tick: (0, 0, 0)) == emission


def test_a_finer_clock_sees_what_a_coarser_one_rounds_away():
    """The same distance at three resolutions: each arrival is the light-travel time rounded up to its tick."""
    metres = 1_234_567.0
    seconds = metres / C
    for ticks, grid in [causal.MILLISECOND_METRE, causal.MICROSECOND_METRE, causal.NANOSECOND_MILLIMETRE]:
        clock = Clock(768, ticks, grid)
        n = clock.light_ticks(clock.grid(metres, 0.0, 0.0))
        assert n == math.ceil(seconds * ticks)
        assert 0.0 <= clock.local_seconds(n) - seconds < clock.tick_seconds
