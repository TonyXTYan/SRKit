"""Constants and small builders shared by the tests."""

from __future__ import annotations

import math
import random

import srkit
from srkit import plan, units

#: 1 g as a rate of change of rapidity per Julian year.
G = units.G0_PER_YEAR
ORIGIN = (0.0, 0.0, 0.0)
X_AXIS = (1.0, 0.0, 0.0)
INF = math.inf
#: Distance to Proxima Centauri in light-years.
D_PROXIMA = 4.2465


def proxima_ship() -> srkit.Worldline:
    """At rest at home, then a 1 g flip-and-burn to Proxima along x from t = 0, then at rest there for ever.

    The ship's clock reads zero at the launch. The leading rest matters: without it the first burn would be
    continued into the past, as a ship that had been braking in from infinity.
    """
    out = plan.line(plan.Start(), (D_PROXIMA, 0.0, 0.0), a=G)
    segments = [srkit.Segment.rest(-1.0, 0.0, ORIGIN), *out.segments, plan.rest(out.end, INF).segment]
    return srkit.Worldline(segments, anchor_t=0.0, anchor_tau=0.0)


def unit(rng: random.Random) -> tuple[float, float, float]:
    """A unit vector in a random direction."""
    a, b = rng.uniform(0.0, math.tau), rng.uniform(-1.5, 1.5)
    return (math.cos(a) * math.cos(b), math.sin(a) * math.cos(b), math.sin(b))


def point(rng: random.Random, size: float = 50.0) -> tuple[float, float, float]:
    return (rng.uniform(-size, size), rng.uniform(-size, size), rng.uniform(-size, size))


def scaled(v, k: float) -> tuple[float, float, float]:
    return (v[0] * k, v[1] * k, v[2] * k)


def segments_of_every_kind(rng: random.Random) -> list[srkit.Segment]:
    """One segment of each kind with parameters in the ranges the game uses, each 30 long at most."""
    seg = srkit.Segment
    t0 = rng.uniform(-100.0, 100.0)
    phi0 = rng.uniform(0.5, 4.0)
    a = rng.uniform(0.01, 3.0)
    r, beta = rng.uniform(0.1, 40.0), rng.uniform(0.01, 0.99)
    burn_alpha = rng.uniform(0.05, 1.0)
    return [
        seg.rest(t0, t0 + 30.0, point(rng)),
        seg.inertial(t0, t0 + 30.0, point(rng), scaled(unit(rng), rng.uniform(0.0, 0.999))),
        seg.hyperbolic(t0, t0 + 30.0, point(rng), unit(rng), a, rng.uniform(0.0, 4.0)),
        # braking, ending before the turnaround
        seg.hyperbolic(t0, t0 + 0.8 * math.sinh(phi0) / a, point(rng), unit(rng), -a, phi0),
        seg.circular(t0, t0 + 30.0, point(rng), r, beta / r, a0=rng.uniform(0.0, math.tau)),
        seg.steered(t0, t0 + 30.0, point(rng), scaled(unit(rng), rng.uniform(0.0, 0.99)), unit(rng), a),
        # a rocket that spends at most half its mass: the burn-out is at proper time w / alpha0
        seg.burn(t0, t0 + 0.5 * 0.7 / burn_alpha, point(rng), unit(rng), rng.uniform(0.0, 1.0), burn_alpha, 0.7),
    ]
