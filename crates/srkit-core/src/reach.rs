//! Where a ship can come to rest, given how much longer it may fly.

use crate::Error;

/// The point where a ship can stop after a given proper time under thrust.
///
/// The ship moves at rapidity `phi` (along `+x`) and can hold proper acceleration `a` for a further proper time
/// `budget`. It burns in a direction fixed aboard at angle `heading` from its velocity, then turns and brakes to
/// rest along its velocity at that moment. The burn time `t1` is the root of
/// `t1 + phi1(t1) / a = budget`, where `phi1` is the rapidity after the burn; the left side increases, so the
/// root is unique. Returns the stopping point as `(along, across)` relative to the ship's position and line of
/// flight, or `None` when the budget cannot even cancel the ship's present speed.
pub fn reach(phi: f64, a: f64, budget: f64, heading: f64) -> Result<Option<[f64; 2]>, Error> {
    if !(phi.is_finite() && phi >= 0.0)
        || !(a.is_finite() && a > 0.0)
        || !(budget.is_finite() && budget >= 0.0)
        || !heading.is_finite()
    {
        return Err(Error::InvalidInput);
    }
    let (g0, gb0) = (libm::cosh(phi), libm::sinh(phi));
    let (c, sn) = (libm::cos(heading), libm::sin(heading));
    let fd = |t1: f64| {
        let (sh, ch) = (libm::sinh(a * t1), libm::cosh(a * t1));
        let ua = g0 * sh * c + gb0 * ch;
        let us = sh * sn;
        let um = libm::hypot(ua, us);
        let dum = if um > 0.0 {
            a * (ua * (g0 * ch * c + gb0 * sh) + us * ch * sn) / um
        } else {
            0.0
        };
        (
            t1 + libm::asinh(um) / a - budget,
            1.0 + dum / libm::hypot(1.0, um) / a,
        )
    };
    if fd(0.0).0 > 1e-12 {
        return Ok(None);
    }
    let (mut lo, mut hi) = (0.0, budget);
    let mut t1 = ((budget - phi / a) / 2.0).max(0.0).min(budget);
    for _ in 0..80 {
        let (f, d) = fd(t1);
        if !f.is_finite() || !d.is_finite() {
            return Err(Error::Overflow);
        }
        if f < 0.0 {
            lo = t1;
        } else {
            hi = t1;
        }
        if f.abs() < 1e-12 * (1.0 + budget) || hi - lo < 1e-13 * (1.0 + budget) {
            break;
        }
        let candidate = if d > 1e-9 { t1 - f / d } else { f64::NAN };
        t1 = if candidate > lo && candidate < hi {
            candidate
        } else {
            (lo + hi) / 2.0
        };
    }
    let (sh, ch) = (libm::sinh(a * t1), libm::cosh(a * t1));
    let k = (g0 - 1.0) * sh * c + gb0 * ch;
    let tm = sh / a;
    let half = libm::sinh(a * t1 / 2.0);
    let hm = 2.0 * half * half / a;
    let (ua, us) = (sh * c + k, sh * sn);
    let um = libm::hypot(ua, us);
    // braking distance (cosh phi1 - 1) / a
    let db = (um / (libm::hypot(1.0, um) + 1.0)) * um / a;
    let out = [
        g0 * hm * c + gb0 * tm + if um > 0.0 { ua * db / um } else { 0.0 },
        hm * sn + if um > 0.0 { us * db / um } else { 0.0 },
    ];
    if out[0].is_finite() && out[1].is_finite() {
        Ok(Some(out))
    } else {
        Err(Error::Overflow)
    }
}

/// Proper time for a rest-to-rest flip-and-burn over distance `d` at proper acceleration `a`, with no coast:
/// `(2 / a) acosh(1 + a d / 2)`.
pub fn flip_burn_tau(d: f64, a: f64) -> f64 {
    if !(d > 0.0) {
        return 0.0;
    }
    2.0 / a * crate::plan::acosh1p(a * d / 2.0)
}

/// The distance a rest-to-rest flip-and-burn covers in proper time `tau` at proper acceleration `a`:
/// `(2 / a)(cosh(a tau / 2) - 1)`. The radius of the "reach from rest" isochrone.
pub fn flip_burn_distance(tau: f64, a: f64) -> f64 {
    2.0 / a * crate::plan::cosh_m1(0.5 * a * tau)
}
