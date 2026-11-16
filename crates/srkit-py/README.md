# srkit for Python

Python bindings for [`srkit-core`](../srkit-core): special relativity in 3+1 Minkowski spacetime, built from
closed forms. Worldlines are piecewise analytic, what an observer sees is a light-cone query, and stored time is
an exact integer.

Status: experimental and private. The licence of SRKit is not decided, so the package declares none and carries
the `Private :: Do Not Upload` classifier. Do not publish wheels.

## Install from source

You need Python 3.9 or later and Rust through [rustup](https://rustup.rs). The Rust toolchain is the one pinned
by `srkit/rust-toolchain.toml` (1.85.1), which rustup picks up from this directory. PyPI and crates.io must be
reachable for the first build.

```sh
cd srkit/crates/srkit-py
python3 -m venv .venv
. .venv/bin/activate
pip install maturin pytest numpy     # NumPy is optional: see "Many times at once"
maturin develop --release            # builds the extension and installs `srkit` into the venv
python -m pytest -q                  # the tests read ../../../tools/physics-checks/oracle.csv
```

`maturin build --release --locked` writes a wheel to `target/wheels/` instead
(`srkit-0.2.0-cp39-abi3-manylinux_2_34_x86_64.whl`, about 0.5 MB, on the Linux machine it was last built on; the
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

Geometry is floating point, which is right for drawing and wrong for deciding. Stored times are whole
milliseconds in 384 or 512 bits and convert to and from Python `int` without loss; `srkit.causal` decides
arrivals in integer arithmetic.

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
- **Exact quantities are `int`.** Milliseconds and metres in `Time384`, `Time512`, the durations and
  `srkit.causal` are Python ints; floats are refused there.

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
| `Time384`, `Time512`, `Duration384`, `Duration512` | Exact integer milliseconds: `+`, `-`, comparisons, `% period_ms`, `local_seconds_since` |
| `causal` | Exact causal decisions: `light_ticks`, `first_arrival`, `arrived`, `classify`, `settle_arrival`, `grid_from_metres` |
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
| `ValueError` | Invalid or out-of-range input: a speed at or above c, a direction that is not a unit vector, a NaN, a run too short for the change of speed asked for, segments that do not join (`check_continuity`) |
| `OverflowError` | A result too large for its representation: a time past `Time384.MAX`, a displacement beyond the signed 128-bit grid |
| `srkit.NoIntersection` | A light cone does not meet a segment within the range given, or a worldline at all (`seen_from`) |
| `srkit.PrecisionLoss` | The answer cannot be given to the promised precision: an exact time difference beyond 2^53 ms asked for as a float. A subclass of `ArithmeticError` |
| `TypeError` | The wrong kind of argument: a float where an exact `int` is needed, times of two different widths |

## Tests

`python -m pytest -q` runs:

- the physics oracle: the 47 rows of `tools/physics-checks/oracle.csv` that are special relativity, optics or
  exact timekeeping (the same rows as the Rust core's oracle test), recomputed through the Python API at each
  row's tolerance;
- identities through the bindings: proper time inverts coordinate time, light-cone roots lie on the cone,
  worldline queries agree with segment queries, vectorised methods equal the scalar ones bit for bit, with and
  without NumPy;
- exact time: Python ints round-trip at both widths, ordering, modular phase at huge epochs, `arrived` flips
  exactly at `first_arrival`;
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
