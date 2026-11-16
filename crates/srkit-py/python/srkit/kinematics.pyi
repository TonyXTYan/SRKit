"""Type stubs for srkit.kinematics: boosts, velocity composition, rapidity, Doppler factors and aberration."""

from typing import Literal, Protocol, SupportsFloat

_Vec3 = tuple[float, float, float]

class _Vec3Like(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> SupportsFloat: ...

def gamma(beta: _Vec3Like) -> float:
    """Lorentz factor of a velocity (a vector, as a fraction of c)."""

def gamma_of_speed(b: float) -> float:
    """Lorentz factor of a speed 0 <= b < 1."""

def rapidity_of_speed(b: float) -> float:
    """Rapidity of a speed: atanh(b)."""

def speed_of_rapidity(phi: float) -> float:
    """Speed of a rapidity, as a fraction of c: tanh(phi)."""

def speed_shortfall(phi: float) -> float:
    """1 - beta for a rapidity, without cancellation."""

def kinetic_per_rest_energy(phi: float) -> float:
    """Kinetic energy per unit rest energy, gamma - 1, from a rapidity."""

def boost(beta: _Vec3Like, dt: float, dx: _Vec3Like) -> tuple[float, _Vec3]:
    """Boost a displacement (dt, dx) into the frame moving at velocity beta. Returns (dt', dx')."""

def boost_rapidity(rapidity: _Vec3Like, dt: float, dx: _Vec3Like) -> tuple[float, _Vec3]:
    """Boost a displacement (dt, dx) by a rapidity vector. Returns (dt', dx')."""

def compose_velocities(v: _Vec3Like, u: _Vec3Like) -> _Vec3:
    """The lab velocity of a body that moves at u in a frame which itself moves at v."""

def relative_velocity(v: _Vec3Like, u: _Vec3Like) -> _Vec3:
    """The velocity of a body moving at u in the lab, as measured in a frame moving at v."""

def doppler_factor(
    emitter_velocity: _Vec3Like, receiver_velocity: _Vec3Like, photon_direction: _Vec3Like
) -> float:
    """Received frequency divided by emitted frequency."""

def doppler_of_rapidity(phi: float) -> float:
    """Head-on Doppler factor for a rapidity: exp(phi)."""

def aberrate(direction: _Vec3Like, beta: _Vec3Like) -> _Vec3:
    """The unit direction a photon travels in, as measured in a frame moving at velocity beta."""

def aberrate_angle(theta: float, phi: float) -> float:
    """The seen angle from the direction of motion of a source at rest-frame angle theta (radians)."""

def interval2(t_a: float, x_a: _Vec3Like, t_b: float, x_b: _Vec3Like) -> float:
    """The invariant interval squared between two events (positive when timelike)."""

def separation(
    t_a: float, x_a: _Vec3Like, t_b: float, x_b: _Vec3Like
) -> Literal["timelike", "null", "spacelike"]:
    """How two events are related, decided in floating point."""

def can_reach(t_a: float, x_a: _Vec3Like, t_b: float, x_b: _Vec3Like) -> bool:
    """Whether light (or anything slower) leaving event a can reach event b."""

def null_coordinates(t: float, x: float) -> tuple[float, float]:
    """Null coordinates along one axis: (u, v) = (t - x, t + x)."""

def from_null_coordinates(u: float, v: float) -> tuple[float, float]:
    """Back from null coordinates: (t, x)."""

__all__ = [
    "aberrate",
    "aberrate_angle",
    "boost",
    "boost_rapidity",
    "can_reach",
    "compose_velocities",
    "doppler_factor",
    "doppler_of_rapidity",
    "from_null_coordinates",
    "gamma",
    "gamma_of_speed",
    "interval2",
    "kinetic_per_rest_energy",
    "null_coordinates",
    "rapidity_of_speed",
    "relative_velocity",
    "separation",
    "speed_of_rapidity",
    "speed_shortfall",
]
