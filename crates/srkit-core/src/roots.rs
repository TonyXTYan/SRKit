//! Bracketed root finding for monotone functions.
//!
//! Light-cone intersections with a timelike worldline are roots of a strictly increasing function, so a bracket
//! plus a safeguarded Newton step always converges, in a bounded number of evaluations and without a time step.

/// Relative width at which [`newton`] stops: a few units in the last place of the root.
pub const NEWTON_STEP_TOLERANCE: f64 = 4e-16;

/// Root of an increasing function on `[lo, hi]`, given `f(lo) <= 0 <= f(hi)`.
///
/// `f` returns the value and the slope. Each step is Newton's where it stays inside the bracket and a bisection
/// otherwise, so the bracket always shrinks. Stops when a step moves the root by less than
/// [`NEWTON_STEP_TOLERANCE`] of its size (or of 1, whichever is larger), when the bracket is twice that, or after
/// 100 evaluations.
pub fn newton<F>(mut f: F, lo: f64, hi: f64, guess: f64) -> f64
where
    F: FnMut(f64) -> (f64, f64),
{
    let (mut a, mut b) = (lo, hi);
    let mut t = guess.max(a).min(b);
    for _ in 0..100 {
        let (value, slope) = f(t);
        if value > 0.0 {
            b = t;
        } else {
            a = t;
        }
        if value == 0.0 {
            return t;
        }
        let mut next = t - value / slope;
        // NaN compares false, so a non-finite step also falls back to the midpoint.
        if !(next > a && next < b) {
            next = 0.5 * (a + b);
        }
        let scale = t.abs().max(1.0);
        if (next - t).abs() <= NEWTON_STEP_TOLERANCE * scale
            || b - a <= 2.0 * NEWTON_STEP_TOLERANCE * scale
        {
            return next;
        }
        t = next;
    }
    t
}

/// Bisection for an increasing function with `f(lo) < 0 <= f(hi)`, returning the smallest representable `t` it
/// reaches with `f(t) >= 0` (`upper = true`) or the largest with `f(t) <= 0` (`upper = false`, for which the
/// bracket is `f(lo) <= 0 < f(hi)`).
///
/// Used where the function can be flat (light running alongside a light-like segment), where Newton's method has
/// no slope to follow.
pub fn bisect<F>(mut f: F, lo: f64, hi: f64, upper: bool) -> f64
where
    F: FnMut(f64) -> f64,
{
    let (mut a, mut b) = (lo, hi);
    for _ in 0..200 {
        if b - a <= NEWTON_STEP_TOLERANCE * b.abs().max(1.0) {
            break;
        }
        let m = 0.5 * (a + b);
        let v = f(m);
        let go_right = if upper { v < 0.0 } else { v <= 0.0 };
        if go_right {
            a = m;
        } else {
            b = m;
        }
    }
    if upper {
        b
    } else {
        a
    }
}

/// The non-negative `u` with `u² = |w0 + v u|²`: how long light from a point takes to meet a body moving at
/// constant velocity `v` (`|v| < 1`), where `w0` is the body's position relative to the point at the moment of
/// emission. Written to avoid cancellation on both branches. Arguments are the dot products `w0·w0`, `w0·v` and
/// `v·v`.
pub fn light_meets_line(ww: f64, wv: f64, vv: f64) -> f64 {
    if ww == 0.0 {
        return 0.0;
    }
    let s = libm::sqrt(wv * wv + (1.0 - vv) * ww);
    if wv >= 0.0 {
        (wv + s) / (1.0 - vv)
    } else {
        ww / (s - wv)
    }
}
