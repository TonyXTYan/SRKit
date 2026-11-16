"""SRKit: special relativity in 3+1 Minkowski spacetime, built from closed forms.

* Worldlines are piecewise analytic. A `Segment` is rest, inertial motion, constant proper acceleration (along
  the motion, or at any angle fixed aboard), a constant-thrust burn that spends mass, a circle at constant
  speed, or a light-like jump. Position, velocity and proper time are evaluated directly at any coordinate
  time; nothing is integrated step by step.
* What an observer sees is a light-cone query. `Worldline.send` finds the moment on a worldline whose light
  reaches an observer, and `Worldline.recv` the moment light from somewhere reaches the worldline.
* Stored time is an exact integer. `Time384` and `Time512` count whole milliseconds and convert to and from
  Python `int` without loss; `srkit.causal` decides "has this signal arrived by tick n?" in integers.

Units
-----
Geometry uses natural units with c = 1: pick any unit of time and measure lengths in the distance light covers
in it (seconds and light-seconds, Julian years and light-years). Velocities are fractions of c, rapidities
have no unit and proper accelerations are rates of change of rapidity per unit of time (1 g is
`srkit.units.G0_PER_YEAR` = 1.0323 per year). `srkit.units` holds the exact SI constants and
`srkit.units.Units` converts to and from SI.

Positions, velocities and directions are `(x, y, z)`: any sequence of three numbers is accepted (a tuple, a
list, a NumPy array) and tuples of floats are returned. An event is given as a time and a position.

Contents
--------
Classes: `Segment`, `Worldline`, `Trip`, `FlipBurn`, `CircularOrbit`, `Time384`, `Time512`, `Duration384`,
`Duration512`, and the small records they return (`TauState`, `Radar`, `RindlerHorizon`, `TripState`, `Leg`,
`OrbitState`).

Functions: `reach`, `flip_burn_tau`, `flip_burn_distance`.

Modules: `units` (constants and conversion), `plan` (itineraries), `kinematics` (boosts, rapidity, Doppler,
aberration), `optics` (what an observer sees; telescopes and beams), `rocket` (the rocket equation and sails),
`kepler` (Kepler's equation), `causal` (exact causal decisions).

Errors
------
Invalid or out-of-range input raises `ValueError`; a result too large for its representation raises
`OverflowError`; a light cone that misses raises `NoIntersection`; an answer that cannot be given to the
promised precision raises `PrecisionLoss`.

NumPy is optional. The vectorised worldline methods (`at_many`, `tau_many`, `send_many`, `recv_many`) return
NumPy arrays when NumPy can be imported at the time of the call and plain lists otherwise.
"""

import sys as _sys
from array import array as _array

from . import _srkit
from ._srkit import (
    CircularOrbit,
    Duration384,
    Duration512,
    FlipBurn,
    Leg,
    NoIntersection,
    OrbitState,
    PrecisionLoss,
    Radar,
    RindlerHorizon,
    Segment,
    TauState,
    Time384,
    Time512,
    Trip,
    TripState,
    flip_burn_distance,
    flip_burn_tau,
    reach,
)

__version__ = _srkit.__version__

# The submodules are built by the extension. Registering them lets `import srkit.units` and
# `from srkit.plan import line` work as for any package.
causal = _srkit.causal
kepler = _srkit.kepler
kinematics = _srkit.kinematics
optics = _srkit.optics
plan = _srkit.plan
rocket = _srkit.rocket
units = _srkit.units
for _module in (causal, kepler, kinematics, optics, plan, rocket, units):
    _sys.modules[_module.__name__] = _module
del _module


def _numpy():
    """NumPy if it can be imported now, otherwise None."""
    try:
        import numpy
    except ImportError:
        return None
    return numpy


def _materialise(values):
    """A one-shot iterator (a generator, say) as a list; anything else unchanged."""
    return list(values) if hasattr(values, "__next__") else values


def _pack_times(ts):
    """Times as a buffer of native doubles: `(data, numpy or None, shape)`."""
    ts = _materialise(ts)
    np = _numpy()
    if np is not None:
        arr = np.asarray(ts, dtype=np.float64)
        # tobytes() copies in C order whatever the array's layout in memory
        return arr.tobytes(), np, arr.shape
    try:
        data = ts if isinstance(ts, _array) and ts.typecode == "d" else _array("d", ts)
    except TypeError:
        raise TypeError("expected a flat sequence of numbers (times)") from None
    return data.tobytes(), None, (len(data),)


def _pack_events(events):
    """Events as a buffer of four native doubles each (t, x, y, z): `(data, numpy or None, shape)`."""
    events = _materialise(events)
    np = _numpy()
    if np is not None:
        arr = np.asarray(events, dtype=np.float64)
        if arr.ndim >= 2 and arr.shape[-1] == 4:
            shape = arr.shape[:-1]
        elif arr.ndim == 1 and arr.size % 4 == 0:
            shape = (arr.size // 4,)
        else:
            raise ValueError(
                "events must have four numbers each (t, x, y, z): an array of shape (n, 4), "
                f"or a flat one of length 4 n; got shape {arr.shape}"
            )
        return arr.tobytes(), np, shape
    if isinstance(events, _array) and events.typecode == "d":
        data = events
    else:
        items = list(events)
        if items and hasattr(items[0], "__len__"):
            flat = []
            for event in items:
                if len(event) != 4:
                    raise ValueError("each event must be four numbers: (t, x, y, z)")
                flat.extend(event)
            items = flat
        try:
            data = _array("d", items)
        except TypeError:
            raise TypeError(
                "expected events as a flat sequence of numbers or a sequence of (t, x, y, z)"
            ) from None
    if len(data) % 4:
        raise ValueError(
            f"events must have four numbers each (t, x, y, z); got {len(data)} numbers"
        )
    return data.tobytes(), None, (len(data) // 4,)


def _unpack(buffer, np, shape):
    """A buffer of native doubles as a NumPy array of the given shape, or as a flat list."""
    if np is not None:
        return np.frombuffer(buffer, dtype=np.float64).reshape(shape)
    out = _array("d")
    out.frombytes(buffer)
    return out.tolist()


class Worldline(_srkit.Worldline):
    __doc__ = _srkit.Worldline.__doc__
    __slots__ = ()

    def at_many(self, ts):
        """Positions and velocities at many coordinate times, as `(xs, vs)`.

        `ts` is any sequence or buffer of numbers (a list, an `array.array`, a NumPy array of any shape).
        With NumPy importable, `xs` and `vs` are float arrays of shape `ts.shape + (3,)`. Without it they
        are flat lists of three floats per time: `[x0, y0, z0, x1, y1, z1, ...]`.

        Each entry equals what `at` returns for that time.
        """
        data, np, shape = _pack_times(ts)
        xs, vs = self._at_many(data)
        return _unpack(xs, np, shape + (3,)), _unpack(vs, np, shape + (3,))

    def tau_many(self, ts):
        """The body's proper time at many coordinate times.

        `ts` is any sequence or buffer of numbers. With NumPy importable the result is a float array of the
        same shape as `ts`; without it, a list of floats.
        """
        data, np, shape = _pack_times(ts)
        return _unpack(self._tau_many(data), np, shape)

    def recv_many(self, events):
        """When light from each of many events reaches this worldline: `recv` for each event.

        `events` holds four numbers per event, `t, x, y, z`: a NumPy array of shape `(n, 4)` (or any shape
        ending in 4), a sequence of `(t, x, y, z)`, or a flat sequence of `4 n` numbers. With NumPy importable
        the result is a float array with one coordinate time per event; without it, a list of floats. An
        entry is `inf` where the light never arrives.
        """
        data, np, shape = _pack_events(events)
        return _unpack(self._cone_many(data, False), np, shape)

    def send_many(self, events):
        """The moments on this worldline seen from each of many events: `send` for each event.

        `events` is as for `recv_many`. The result holds one coordinate time per event (a NumPy float array
        if NumPy is importable, otherwise a list of floats). An entry is `-inf` where no such moment exists.
        """
        data, np, shape = _pack_events(events)
        return _unpack(self._cone_many(data, True), np, shape)


__all__ = [
    "CircularOrbit",
    "Duration384",
    "Duration512",
    "FlipBurn",
    "Leg",
    "NoIntersection",
    "OrbitState",
    "PrecisionLoss",
    "Radar",
    "RindlerHorizon",
    "Segment",
    "TauState",
    "Time384",
    "Time512",
    "Trip",
    "TripState",
    "Worldline",
    "__version__",
    "causal",
    "flip_burn_distance",
    "flip_burn_tau",
    "kepler",
    "kinematics",
    "optics",
    "plan",
    "reach",
    "rocket",
    "units",
]
