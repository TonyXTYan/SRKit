# srkit for Python

Python bindings for [`srkit-core`](../srkit-core): special relativity in 3+1 Minkowski spacetime, built from
closed forms. Worldlines are piecewise analytic, what an observer sees is a light-cone query, and stored time is
an exact integer.

Status: experimental. MIT licence, as the rest of SRKit. The published name is not settled, so the package
carries the `Private :: Do Not Upload` classifier: do not publish wheels yet.

## Install from source

You need Python 3.9 or later and Rust through [rustup](https://rustup.rs). The Rust toolchain is the one pinned
by `rust-toolchain.toml` (1.85.1), which rustup picks up from this directory. PyPI and crates.io must be
reachable for the first build.

```sh
cd crates/srkit-py
python3 -m venv .venv
. .venv/bin/activate
pip install maturin pytest numpy     # NumPy is optional: see "Many times at once"
maturin develop --release            # builds the extension and installs `srkit` into the venv
python -m pytest -q                  # the tests read ../../oracle/oracle.csv
```

`maturin build --release --locked` writes a wheel to `target/wheels/` instead
(`srkit-0.3.0-cp39-abi3-manylinux_2_34_x86_64.whl`, about 0.5 MB, on the Linux machine it was last built on; the
`manylinux` tag follows the build host's C library). It is an `abi3` wheel: one file serves CPython 3.9 and
every later version.

The crate is its own cargo workspace (the `srkit` workspace excludes it, because it needs a Python toolchain)
and depends on the core by path, so build it from a checkout of the repository. Build the extension with
maturin: a plain `cargo build` links the Python library in, which an extension module must not do.

## Example: a 1 g trip to Proxima, seen from home

```python
import math
from srkit import Segment, Worldline, plan, units

g = units.G0_PER_YEAR                          # 1 g as a rate of change of rapidity: 1.0323 per year
home, proxima = (0, 0, 0), (4.2465, 0, 0)      # light-years

# Plan the trip, rest to rest: boost to the midpoint, turn, brake.
trip = plan.line(plan.Start(t=0.0, tau=0.0, x=home), proxima, a=g)
print(f"{trip.end.t:.3f} yr at home, {trip.end.tau:.3f} yr aboard, peak {math.tanh(trip.phi_peak):.4f} c")

# Its worldline: waiting at home, the two burns, then parked at Proxima for ever.
ship = Worldline(
    [Segment.rest(-1.0, 0.0, home), *trip.segments, plan.rest(trip.end, math.inf).segment],
    anchor_t=0.0,                              # the ship's clock reads zero at the launch
)
ship.check_continuity()

# What does an observer at home see three years after the launch? The light arriving then left the ship earlier.
t_seen, x, v = ship.seen_from(3.0, home)
print(f"home at 3 yr sees the ship as it was at {t_seen:.3f} yr: "
      f"{x[0]:.3f} ly out at {v[0]:.3f} c, its clock reading {ship.tau(t_seen):.3f} yr")
```

```text
5.873 yr at home, 3.542 yr aboard, peak 0.9497 c
home at 3 yr sees the ship as it was at 1.866 yr: 1.134 ly out at 0.888 c, its clock reading 1.366 yr
```

The leading `Segment.rest` matters. A worldline's first segment is continued into the infinite past, so without
it the first burn would describe a ship that had been braking in from infinity, not one waiting at home.

## Example: an exact arrival a googol of milliseconds from the epoch

Geometry is floating point, which is right for drawing and wrong for deciding. Stored times are whole ticks in
384, 512, 768 or 1024 bits (`Time384`, `Time512`, `Time768`, `Time1024`) and convert to and from Python `int`
without loss; `srkit.causal` decides arrivals in integer arithmetic. A tick is a millisecond and the grid a
metre unless you choose otherwise (next section).

```python
from srkit import Time512, causal, units

launch = Time512(10**100)                                        # 10^100 ms after the epoch
dx = causal.grid_from_metres(4.2465 * units.LIGHT_YEAR, 0.0, 0.0)   # whole metres, as ints
arrival = causal.first_arrival(launch, dx)                       # first tick at which the signal is there
print(arrival - launch)                                          # the light-travel time in ms, rounded up
assert causal.arrived(launch, arrival, dx) and not causal.arrived(launch, arrival - 1, dx)
assert int(arrival) == 10**100 + 134_009_348_401                 # exact: no float was involved
assert causal.light_ticks((299_792_458, 0, 0)) == 1000           # one light-second is exactly 1000 ticks
assert causal.light_ticks((299_792_459, 0, 0)) == 1001           # and one metre more is 1001
```

```text
134009348401
```

## Example: choosing the width, the tick and the grid

A time counts ticks and does not know how long a tick is. `causal.Clock(bits, ticks_per_second, grid_per_metre)`
holds the three choices together: the width of the stored integer (384, 512, 768 or 1024 bits), the tick (from
a second down to a picosecond) and the grid spacing (from a metre down to a micrometre). Its methods are the
same exact checks at that resolution. Displacements are three ints of grid steps.

```python
from srkit import causal

clock = causal.Clock(768, *causal.NANOSECOND_MILLIMETRE)     # 768 bits, nanosecond ticks, millimetre grid
print(clock, f"reaches 10^{clock.reach_log10_years:.1f} years")
sent = clock.time(3 * 10**39)                                # a Time768: 3e39 ns after the epoch
away = clock.grid(299_792_458.0, 0.0, 0.0)                   # one light-second, in whole millimetres
heard = clock.first_arrival(sent, away)                      # first tick at which the signal is there
print(heard - sent)                                          # a second of nanoseconds, exactly
assert clock.arrived(sent, heard, away) and not clock.arrived(sent, heard - 1, away)
assert clock.classify(heard - sent, away) == "null"          # exactly on the light cone
assert heard - sent == clock.seconds(1) and clock.local_seconds(heard - sent) == 1.0

default = causal.Clock()                                     # 384 bits, millisecond ticks, metre grid
assert default.first_arrival(default.time(0), (299_792_458, 0, 0)) == default.time(1000)
```

```text
Clock(bits=768, ticks_per_second=1000000000, grid_per_metre=1000) reaches 10^214.7 years
1000000000
```

The defaults are `Clock(384, 1000, 1)`: milliseconds on a metre grid. The module-level functions
(`causal.first_arrival`, `causal.arrived`, ...) are that default at the width of the times they are given, so
`causal.first_arrival(t, dx)` is `Clock(t.BITS).first_arrival(t, dx)`. `MILLISECOND_METRE`, `MICROSECOND_METRE`
and `NANOSECOND_MILLIMETRE` are `(ticks_per_second, grid_per_metre)` pairs to unpack into `Clock`.

Three things to know when choosing:

- **The width and the tick set the reach.** `clock.reach_log10_years` is the longest time the clock holds, as
  a power of ten in years: 105.1 for 384 bits of milliseconds, 99.1 for 384 bits of nanoseconds (just short of
  a googol), 214.7 for 768 bits of nanoseconds.
- **The grid must be much finer than the distance light covers in a tick**, or rounding a position to the grid
  moves an arrival by whole ticks. `clock.light_steps_per_tick` should be in the hundreds or more: a
  nanosecond tick wants a millimetre grid (300 steps), not a metre one (0.3).
- **A time does not carry its clock.** `Time768(5)` is five ticks, whatever a tick is. Keep one clock per
  application and use it for every call. The methods on the time classes that give seconds
  (`local_seconds_since`, `phase_seconds`, `stationary_arrival`, ...) read a tick as a millisecond; for other
  ticks use the clock's methods (`local_seconds_between`, `phase_seconds`, `first_arrival`, ...). Times of
  different widths never mix: that is a `TypeError`.

## Units and conventions

- **Natural units, c = 1.** Pick a unit of time and measure lengths in the distance light covers in it: years
  and light-years, seconds and light-seconds. Velocities are fractions of c, rapidities have no unit, and a
  proper acceleration is a rate of change of rapidity per unit of time. `srkit.units` holds the exact SI
  constants (`C`, `JULIAN_YEAR`, `LIGHT_YEAR`, `AU`, `PARSEC`, `G0`, `G0_PER_YEAR`, ...) and
  `units.Units` converts: `Units.YEARS.accel_from_si(9.80665)` is 1.0323.
- **Vectors** are `(x, y, z)`. Any sequence of three numbers is accepted (a tuple, a list, a NumPy array);
  tuples of floats come back. An event is passed as a time and a position: `recv(event_t, event_x)`.
- **Numbers pass through unchanged.** The bindings check arguments (a velocity below c, a unit vector where one
  is needed, no NaN) but never normalise or round them, so a call gives the same bits as the same call on the
  Rust core.
- **Exact quantities are `int`.** Ticks and grid steps (milliseconds and metres by default) in the time
  classes, the durations and `srkit.causal` are Python ints; floats are refused there.

## What is in the package

| Name | What it is |
|---|---|
| `Segment` | One closed-form piece of a worldline. Constructors: `rest`, `inertial`, `hyperbolic`, `steered`, `burn`, `circular`, `null`, `jump`. Methods: `at(t)`, `dtau(t)`, `t_of_tau(s)`, `at_tau(s)`, `vmax()`, `recv(...)`, `send(...)`, `horizon_margin(...)`, `rindler_horizon()` |
| `Worldline` | Segments in order, with proper time carried along the chain: `at`, `tau`, `t_of_tau`, `at_tau`, `recv`, `send`, `seen_from`, `radar`, `light_cone_coordinates`, `check_continuity`, and the vectorised `at_many`, `tau_many`, `recv_many`, `send_many` |
| `plan` | Itineraries that join continuously: `line`, `brake`, `orbit_entry`, `orbit_exit`, `plane_through`, and single pieces (`rest`, `hyp_piece`, `lin_piece`, `null_piece`, `circ_piece`) |
| `Trip`, `FlipBurn`, `reach` | Rest-to-rest rides in null coordinates, exact at any Lorentz factor; where a ship can stop within a proper-time budget |
| `CircularOrbit` | A circle flown under thrust: period, clock rates, Sagnac gap, Thomas precession |
| `kinematics`, `optics` | Boosts, rapidity, velocity composition, Doppler factors, aberration; what a moving observer sees; telescope and beam formulas |
| `rocket`, `kepler` | The relativistic rocket equation and light sails; Kepler's equation |
| `Time384`, `Time512`, `Time768`, `Time1024` and `Duration384`, `Duration512`, `Duration768`, `Duration1024` | Exact integer ticks, a millisecond each by default: `+`, `-`, comparisons, `% period`, `local_seconds_since` |
| `causal` | Exact causal decisions with millisecond ticks on a metre grid: `light_ticks`, `first_arrival`, `arrived`, `classify`, `settle_arrival`, `grid_from_metres` |
| `causal.Clock` | A width, a tick and a grid spacing: `time`, `duration`, `seconds`, `grid`, the same decisions at that resolution, `local_seconds`, `local_seconds_between`, `phase_seconds`, `reach_log10_years` |
| `units` | Exact SI constants and `Units` |

Every public name has a docstring (`help(srkit)`), and the package ships type stubs (`py.typed`).

## Light cones

`Worldline.send(event_t, event_x)` is the time on the worldline whose light reaches the event: what an observer
there sees of the body. `Worldline.recv(event_t, event_x)` is the time at which light from the event reaches
the worldline. Both are unique for anything slower than light. On a worldline they return `-inf` and `inf` when
there is no such moment (a ship accelerating away for ever outruns light from behind its horizon). On a single
`Segment` they search a range, `lo` to `hi`, and raise `srkit.NoIntersection` if the light does not arrive
within it.

## Many times at once

The vectorised methods take any sequence or buffer of numbers and call into Rust once:

```python
import numpy as np
from srkit import Segment, Worldline

at_rest = Segment.rest(-1.0, 0.0, (0, 0, 0))
burning = Segment.hyperbolic(0.0, np.inf, (0, 0, 0), (1, 0, 0), alpha=1.0)
ship = Worldline([at_rest, burning], anchor_t=0.0)
ts = np.linspace(0.0, 5.0, 6)
xs, vs = ship.at_many(ts)                      # two arrays of shape (6, 3)
taus = ship.tau_many(ts)                       # shape (6,): asinh(t) for this ship
seen = ship.send_many([(t, 0.0, 0.0, 0.0) for t in ts])   # events are (t, x, y, z)
assert xs.shape == (6, 3) and np.allclose(taus, np.arcsinh(ts)) and seen.shape == (6,)
```

NumPy is not required. The methods return NumPy arrays when NumPy can be imported at the time of the call and
flat lists of floats otherwise (`at_many` then gives three floats per time).

## Errors

| Raised | When |
|---|---|
| `ValueError` | Invalid or out-of-range input: a speed at or above c, a direction that is not a unit vector, a NaN, a run too short for the change of speed asked for, segments that do not join (`check_continuity`), a clock width or resolution that is not offered |
| `OverflowError` | A result too large for its representation: a time past `Time384.MAX`, a displacement beyond the signed 128-bit grid |
| `srkit.NoIntersection` | A light cone does not meet a segment within the range given, or a worldline at all (`seen_from`) |
| `srkit.PrecisionLoss` | The answer cannot be given to the promised precision: an exact time difference beyond 2^53 ticks asked for as a float. A subclass of `ArithmeticError` |
| `TypeError` | The wrong kind of argument: a float where an exact `int` is needed, times of two different widths, a time that is not of a clock's width |

## Tests

`python -m pytest -q` runs:

- the physics oracle: the 47 rows of `oracle/oracle.csv` that are special relativity, optics or
  exact timekeeping (the same rows as the Rust core's oracle test), recomputed through the Python API at each
  row's tolerance;
- identities through the bindings: proper time inverts coordinate time, light-cone roots lie on the cone,
  worldline queries agree with segment queries, vectorised methods equal the scalar ones bit for bit, with and
  without NumPy;
- exact time: Python ints round-trip at all four widths, ordering, modular phase at huge epochs, `arrived`
  flips exactly at `first_arrival`;
- clocks: one light-second is exactly one second of ticks at seven resolutions and four widths, the widths
  agree with each other and with Python's own integers, and the module-level functions equal the default
  clock;
- error mapping, and that the type stubs list exactly the names, members and signatures the module has
  (with `mypy.stubtest` as a second opinion where mypy is installed);
- the examples in this file, including what they print.

Tests that need NumPy are skipped where it is not installed; the rest check the list-returning paths.

Rust checks, from this directory: `cargo clippy --all-targets -- -D warnings` and `cargo fmt -- --check`.

## Not here yet

- Astropy `Quantity` adapters (named in the package plan as optional).
- `srkit-flrw` and `srkit-gr` (cosmology and closed-form gravity) are separate crates and are not bound.
- From the core: `frames.RestFrame` and four-vector arithmetic (use `kinematics.boost` and
  `kinematics.boost_rapidity`), and the root finders in `roots`.
- Wheels for other platforms: only a local Linux build has been made and tested.
