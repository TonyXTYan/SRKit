"""Type stubs for srkit.causal: exact causal decisions on stored integer times.

Times are whole milliseconds (Time384, Time512); displacements are three ints of whole metres.
"""

from typing import Callable, Literal, Protocol, SupportsIndex, TypeVar

from . import Duration384, Duration512, Time384, Time512

class _GridLike(Protocol):
    """A displacement: any sequence of three ints, in metres."""

    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> SupportsIndex: ...

_Time = TypeVar("_Time", Time384, Time512)

C_METRES_PER_SECOND: int
"""Speed of light in metres per second, exactly."""

def light_ticks(dx: _GridLike) -> int:
    """The light-travel time over dx, rounded up to whole milliseconds. Exact."""

def first_arrival(emission: _Time, dx: _GridLike) -> _Time:
    """The first tick at which a signal emitted at `emission` has reached a fixed point dx away."""

def arrived(emission: _Time, now: _Time, dx: _GridLike) -> bool:
    """Whether a signal emitted at `emission` has reached a point dx away by the tick `now`. Exact."""

def classify(
    dt: int | Duration384 | Duration512, dx: _GridLike
) -> Literal["timelike", "null", "spacelike"]:
    """Classify the separation between two events dt milliseconds and dx metres apart. Exact."""

def grid_from_metres(x: float, y: float, z: float) -> tuple[int, int, int]:
    """Round a displacement in metres, given as floats, to the one-metre grid."""

def settle_arrival(
    emission: _Time,
    estimate: _Time,
    displacement_at: Callable[[_Time], _GridLike],
    guard: int = 8,
) -> _Time:
    """The first tick at which a signal reaches a moving receiver, settled exactly from an estimate."""

__all__ = [
    "C_METRES_PER_SECOND",
    "arrived",
    "classify",
    "first_arrival",
    "grid_from_metres",
    "light_ticks",
    "settle_arrival",
]
