"""Type stubs for srkit.plan: itineraries, as sequences of segments that join continuously."""

from typing import Literal, Protocol, SupportsFloat, final

from . import CircularOrbit, Segment

_Vec3 = tuple[float, float, float]

class _Vec3Like(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> SupportsFloat: ...

_Phase = Literal["rest", "boost", "coast", "brake", "orbit", "jump"]

@final
class Start:
    """Where and when a plan (or its next piece) begins."""

    def __new__(cls, t: float = 0.0, tau: float = 0.0, x: _Vec3Like | None = None) -> Start: ...
    @property
    def t(self) -> float: ...
    @property
    def tau(self) -> float: ...
    @property
    def x(self) -> _Vec3: ...

@final
class End:
    """The state at the end of a piece or a plan."""

    @property
    def t(self) -> float: ...
    @property
    def tau(self) -> float: ...
    @property
    def x(self) -> _Vec3: ...
    @property
    def v(self) -> _Vec3: ...
    @property
    def phi(self) -> float: ...
    @property
    def acc(self) -> float: ...
    def start(self) -> Start: ...

@final
class Piece:
    """One piece of a plan: a segment, the proper time at its start, its proper duration and its purpose."""

    @property
    def segment(self) -> Segment: ...
    @property
    def tau0(self) -> float: ...
    @property
    def dtau(self) -> float: ...
    @property
    def phase(self) -> _Phase: ...
    def end(self) -> End: ...

@final
class Line:
    """A straight run and its summary."""

    @property
    def pieces(self) -> tuple[Piece, ...]: ...
    @property
    def segments(self) -> tuple[Segment, ...]: ...
    @property
    def end(self) -> End: ...
    @property
    def n(self) -> _Vec3: ...
    @property
    def dphi(self) -> float: ...
    @property
    def length(self) -> float: ...
    @property
    def phi_peak(self) -> float: ...
    @property
    def coast_length(self) -> float: ...

@final
class Brake:
    """A braking burn to rest and its summary."""

    @property
    def pieces(self) -> tuple[Piece, ...]: ...
    @property
    def segments(self) -> tuple[Segment, ...]: ...
    @property
    def end(self) -> End: ...
    @property
    def dphi(self) -> float: ...

@final
class Plane:
    """An orthonormal pair (e1, e2) and their normal."""

    @property
    def e1(self) -> _Vec3: ...
    @property
    def e2(self) -> _Vec3: ...
    @property
    def normal(self) -> _Vec3: ...

@final
class OrbitElements:
    """The elements of a circle to fly: centre, radius, angular velocity, plane and the angle at the join."""

    def __new__(
        cls,
        centre: _Vec3Like,
        r: float,
        omega: float,
        e1: _Vec3Like | None = None,
        e2: _Vec3Like | None = None,
        a0: float = 0.0,
    ) -> OrbitElements: ...
    @property
    def centre(self) -> _Vec3: ...
    @property
    def r(self) -> float: ...
    @property
    def omega(self) -> float: ...
    @property
    def e1(self) -> _Vec3: ...
    @property
    def e2(self) -> _Vec3: ...
    @property
    def a0(self) -> float: ...
    @property
    def orbit(self) -> CircularOrbit: ...

@final
class OrbitEntry:
    """The pieces that bring a ship onto a circular orbit, and the circle's elements at the join."""

    @property
    def pieces(self) -> tuple[Piece, ...]: ...
    @property
    def segments(self) -> tuple[Segment, ...]: ...
    @property
    def end(self) -> End: ...
    @property
    def orbit(self) -> OrbitElements: ...
    @property
    def dphi(self) -> float: ...
    @property
    def from_inside(self) -> bool: ...

def line(
    start: Start | End,
    x1: _Vec3Like,
    a: float,
    phi_cap: float = ...,
    phi0: float = 0.0,
    phi1: float = 0.0,
    a_brake: float | None = None,
    n_default: _Vec3Like | None = None,
) -> Line:
    """A straight run to x1: boost at a, coast if the rapidity reaches phi_cap, brake at a_brake."""

def brake(start: Start | End, n: _Vec3Like, phi0: float, a: float) -> Brake:
    """Brake from rapidity phi0 along the unit vector n to rest, at proper acceleration a."""

def orbit_exit(start: Start | End, v: _Vec3Like, a: float) -> Brake:
    """Leave a circle while moving at velocity v: brake along the tangent to rest."""

def plane_through(centre: _Vec3Like, x: _Vec3Like, e1_default: _Vec3Like | None = None) -> Plane:
    """The plane through centre that holds the point x and is tilted least from the xy plane."""

def orbit_entry(
    start: Start | End,
    centre: _Vec3Like,
    r: float,
    beta: float,
    a: float,
    phi_cap: float = ...,
    *,
    e1: _Vec3Like | None = None,
    a0: float | None = None,
    a_tangent: float | None = None,
) -> OrbitEntry:
    """Enter a circular orbit of radius r at speed beta about centre, from rest at start."""

def rest(start: Start | End, dtau: float) -> Piece:
    """Stay at rest for proper time dtau."""

def hyp_piece(
    start: Start | End,
    n: _Vec3Like,
    alpha: float,
    phi0: float,
    phi1: float,
    phase: _Phase = "boost",
) -> Piece:
    """Accelerate at alpha along the unit vector n, from rapidity phi0 to phi1."""

def lin_piece(start: Start | End, v: _Vec3Like, dtau: float) -> Piece:
    """Coast at velocity v for proper time dtau."""

def null_piece(start: Start | End, x1: _Vec3Like) -> Piece:
    """A light-like jump to the point x1."""

def circ_piece(start: Start | End, orbit: OrbitElements, dtau: float = ...) -> Piece:
    """Fly the circle for proper time dtau, by default for ever."""

__all__ = [
    "Brake",
    "End",
    "Line",
    "OrbitElements",
    "OrbitEntry",
    "Piece",
    "Plane",
    "Start",
    "brake",
    "circ_piece",
    "hyp_piece",
    "lin_piece",
    "line",
    "null_piece",
    "orbit_entry",
    "orbit_exit",
    "plane_through",
    "rest",
]
