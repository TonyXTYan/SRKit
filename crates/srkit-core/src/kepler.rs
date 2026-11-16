//! Kepler's equation, for bound two-body orbits drawn or used as landmarks.

/// Eccentric anomaly `E` solving `E - e sin E = M` for mean anomaly `M` (radians) and eccentricity `0 <= e < 1`.
///
/// Newton's method from `E = M + e sin M`, at most 16 steps, stopping when a step is below 1e-15. The iteration
/// count is fixed by the inputs, so the result is deterministic.
pub fn eccentric_anomaly(mean_anomaly: f64, e: f64) -> f64 {
    let m = mean_anomaly;
    let mut big_e = m + e * libm::sin(m);
    for _ in 0..16 {
        let d = (big_e - e * libm::sin(big_e) - m) / (1.0 - e * libm::cos(big_e));
        big_e -= d;
        if d.abs() < 1e-15 {
            break;
        }
    }
    big_e
}

/// Position in the orbital plane, in units of the semi-major axis, with the focus at the origin and periapsis
/// along `+x`: `(cos E - e, sqrt(1 - e²) sin E)`.
pub fn position_in_plane(eccentric_anomaly: f64, e: f64) -> (f64, f64) {
    (
        libm::cos(eccentric_anomaly) - e,
        libm::sqrt(1.0 - e * e) * libm::sin(eccentric_anomaly),
    )
}

/// True anomaly from eccentric anomaly.
pub fn true_anomaly(eccentric_anomaly: f64, e: f64) -> f64 {
    2.0 * libm::atan2(
        libm::sqrt(1.0 + e) * libm::sin(0.5 * eccentric_anomaly),
        libm::sqrt(1.0 - e) * libm::cos(0.5 * eccentric_anomaly),
    )
}
