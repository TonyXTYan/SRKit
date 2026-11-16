"""Type stubs for srkit.causal: exact causal decisions on stored integer times.

Times are whole ticks (Time384, Time512, Time768, Time1024); displacements are three ints of whole grid steps.
The module-level functions use a millisecond tick on a metre grid; a Clock holds any width, tick and grid.
"""

from typing import Callable, Literal, Protocol, SupportsIndex, TypeVar, final

from . import (
    Duration384,
    Duration512,
    Duration768,
    Duration1024,
    Time384,
    Time512,
    Time768,
    Time1024,
)

class _GridLike(Protocol):
    """A displacement: any sequence of three ints, in grid steps (metres by default)."""

    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> SupportsIndex: ...

_Time = TypeVar("_Time", Time384, Time512, Time768, Time1024)
_AnyTime = Time384 | Time512 | Time768 | Time1024
_AnyDuration = Duration384 | Duration512 | Duration768 | Duration1024
_Causal = Literal["timelike", "null", "spacelike"]

C_METRES_PER_SECOND: int
"""Speed of light in metres per second, exactly."""

MILLISECOND_METRE: tuple[int, int]
"""(ticks_per_second, grid_per_metre) for a millisecond tick on a metre grid, the default: (1000, 1)."""

MICROSECOND_METRE: tuple[int, int]
"""(ticks_per_second, grid_per_metre) for a microsecond tick on a metre grid: (10**6, 1)."""

NANOSECOND_MILLIMETRE: tuple[int, int]
"""(ticks_per_second, grid_per_metre) for a nanosecond tick on a millimetre grid: (10**9, 1000)."""

@final
class Clock:
    """A stored clock: a width (384, 512, 768 or 1024 bits), a tick and a grid spacing, chosen together.

    Its times are the time class of its width; times of another width are refused with TypeError.
    Displacements are three ints of grid steps.
    """

    def __new__(cls, bits: int = 384, ticks_per_second: int = 1000, grid_per_metre: int = 1) -> Clock: ...
    @property
    def bits(self) -> int: ...
    @property
    def ticks_per_second(self) -> int: ...
    @property
    def grid_per_metre(self) -> int: ...
    @property
    def light_steps_per_tick(self) -> float: ...
    @property
    def tick_seconds(self) -> float: ...
    @property
    def grid_metres(self) -> float: ...
    @property
    def reach_log10_years(self) -> float: ...
    @property
    def Time(self) -> type[Time384] | type[Time512] | type[Time768] | type[Time1024]: ...
    @property
    def Duration(self) -> type[Duration384] | type[Duration512] | type[Duration768] | type[Duration1024]: ...
    def time(self, ticks: int) -> _AnyTime: ...
    def duration(self, ticks: int) -> _AnyDuration: ...
    def seconds(self, n: int) -> _AnyDuration: ...
    def grid(self, x_m: float, y_m: float, z_m: float) -> tuple[int, int, int]: ...
    def light_ticks(self, dx: _GridLike) -> int: ...
    def first_arrival(self, emission: _Time, dx: _GridLike) -> _Time: ...
    def arrived(self, emission: _Time, now: _Time, dx: _GridLike) -> bool: ...
    def classify(self, dt: int | _AnyDuration, dx: _GridLike) -> _Causal: ...
    def settle_arrival(
        self,
        emission: _Time,
        estimate: _Time,
        displacement_at: Callable[[_Time], _GridLike],
        guard: int = 8,
    ) -> _Time: ...
    def local_seconds(self, duration: int | _AnyDuration) -> float: ...
    def local_seconds_between(self, later: _Time, earlier: _Time) -> float: ...
    def phase_seconds(self, t: _AnyTime, period_ticks: int) -> float: ...
    def __hash__(self) -> int: ...
    def __eq__(self, other: object, /) -> bool: ...
    def __ne__(self, other: object, /) -> bool: ...

def light_ticks(dx: _GridLike) -> int:
    """The light-travel time over dx metres, rounded up to whole milliseconds. Exact."""

def first_arrival(emission: _Time, dx: _GridLike) -> _Time:
    """The first millisecond tick at which a signal emitted at `emission` has reached a fixed point dx away."""

def arrived(emission: _Time, now: _Time, dx: _GridLike) -> bool:
    """Whether a signal emitted at `emission` has reached a point dx away by the tick `now`. Exact."""

def classify(dt: int | _AnyDuration, dx: _GridLike) -> _Causal:
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
    "Clock",
    "MICROSECOND_METRE",
    "MILLISECOND_METRE",
    "NANOSECOND_MILLIMETRE",
    "arrived",
    "classify",
    "first_arrival",
    "grid_from_metres",
    "light_ticks",
    "settle_arrival",
]
