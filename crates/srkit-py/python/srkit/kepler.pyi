"""Type stubs for srkit.kepler: Kepler's equation for bound orbits. Angles are in radians."""

def eccentric_anomaly(mean_anomaly: float, e: float) -> float:
    """Eccentric anomaly E solving E - e sin E = M, for eccentricity 0 <= e < 1."""

def position_in_plane(eccentric_anomaly: float, e: float) -> tuple[float, float]:
    """Position in the orbital plane in units of the semi-major axis, focus at the origin."""

def true_anomaly(eccentric_anomaly: float, e: float) -> float:
    """True anomaly from eccentric anomaly."""

__all__ = ["eccentric_anomaly", "position_in_plane", "true_anomaly"]
