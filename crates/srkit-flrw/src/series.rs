//! The Gauss hypergeometric functions the closed forms need.
//!
//! Every integral in this crate reduces to `2F1(a, b; c; -x)` with `x = (Ω_Λ/Ω_m) a³ ≥ 0`. The Pfaff
//! transformation
//!
//! ```text
//! 2F1(a, b; c; z) = (1 - z)^(-b) 2F1(c - a, b; c; z / (z - 1))
//! ```
//!
//! turns each of them into `2F1(1, 1/2; c; y)` with `y = x/(1 + x)` (Abramowitz & Stegun 15.3.4; DLMF 15.8.1).
//! For `x ≤ 1` that argument is at most 1/2. For `x > 1` the crate uses the complementary integral, which has
//! the same shape with `x` replaced by `1/x`, so `y = 1/(1 + x)` is again at most 1/2. The series is therefore
//! only ever summed for `0 ≤ y ≤ 1/2`:
//!
//! ```text
//! 2F1(1, 1/2; c; y) = Σ t_n y^n,   t_0 = 1,   t_(n+1) / t_n = (n + 1/2) / (n + c)
//! ```
//!
//! All terms are positive and `t_n ≤ 1` for the four values of `c` used (7/6, 4/3, 11/6, 2/3), so the tail
//! after degree `N` is below `2 y^(N+1)`. The sum is evaluated by Horner's rule from the highest term down,
//! with a fixed number of terms chosen from `y` alone, so the smallest terms are added first and the result is
//! good to a few units in the last place.

/// Highest degree ever summed: at `y = 1/2` the tail is below `2^-60`.
const TERMS: usize = 60;

/// Term ratios `(6n + 3) / (6n + d)` for `n = 0..TERMS`, which is `(n + 1/2) / (n + d/6)`.
const fn ratios(d: f64) -> [f64; TERMS] {
    let mut out = [0.0; TERMS];
    let mut n = 0;
    while n < TERMS {
        let k = 6.0 * n as f64;
        out[n] = (k + 3.0) / (k + d);
        n += 1;
    }
    out
}

static RATIOS_7_6: [f64; TERMS] = ratios(7.0);
static RATIOS_4_3: [f64; TERMS] = ratios(8.0);
static RATIOS_11_6: [f64; TERMS] = ratios(11.0);
static RATIOS_2_3: [f64; TERMS] = ratios(4.0);

/// Degree that leaves a relative tail below `2^-60` for this `y`.
fn degree(y: f64) -> usize {
    if y <= 0.015_625 {
        10
    } else if y <= 0.125 {
        20
    } else if y <= 0.25 {
        30
    } else {
        TERMS
    }
}

/// `Σ t_n y^n` by Horner's rule, for `0 ≤ y ≤ 1/2`.
fn sum(ratios: &[f64; TERMS], y: f64) -> f64 {
    let mut s = 1.0;
    let mut n = degree(y);
    while n > 0 {
        n -= 1;
        s = 1.0 + y * ratios[n] * s;
    }
    s
}

/// `2F1(1, 1/2; 7/6; y)` for `0 ≤ y ≤ 1/2`: conformal time before matter-Λ equality.
pub(crate) fn f_7_6(y: f64) -> f64 {
    sum(&RATIOS_7_6, y)
}

/// `2F1(1, 1/2; 4/3; y)` for `0 ≤ y ≤ 1/2`: conformal time left after matter-Λ equality.
pub(crate) fn f_4_3(y: f64) -> f64 {
    sum(&RATIOS_4_3, y)
}

/// `2F1(1, 1/2; 11/6; y)` for `0 ≤ y ≤ 1/2`: `∫ a dt` before matter-Λ equality.
pub(crate) fn f_11_6(y: f64) -> f64 {
    sum(&RATIOS_11_6, y)
}

/// `2F1(1, 1/2; 2/3; y)` for `0 ≤ y ≤ 1/2`: `∫ a dt` after matter-Λ equality.
pub(crate) fn f_2_3(y: f64) -> f64 {
    sum(&RATIOS_2_3, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values from mpmath `hyp2f1(1, 1/2, c, y)` at 40 digits, rounded to binary64.
    const REFERENCE: [(f64, [f64; 4]); 6] = [
        (
            0.5,
            [
                1.3397976221021002,
                1.286374656332034,
                1.1912787652628032,
                1.696086648107738,
            ],
        ),
        (
            0.37,
            [
                1.216522994482272,
                1.1849121566413496,
                1.1272452408683205,
                1.4195925520679185,
            ],
        ),
        (
            0.25,
            [
                1.130331066321234,
                1.1123663925084712,
                1.0790079870542748,
                1.2426861517846624,
            ],
        ),
        (
            0.125,
            [
                1.058718841508521,
                1.0510347451975754,
                1.0365492724450034,
                1.1057005414575676,
            ],
        ),
        (
            0.015625,
            [
                1.0067697713420316,
                1.0059189288607473,
                1.0042969772086747,
                1.0118859954648989,
            ],
        ),
        (
            0.001,
            [
                1.0004288683663045,
                1.0003752412523783,
                1.0002728717519864,
                1.0007506756334172,
            ],
        ),
    ];

    #[test]
    fn series_match_mpmath_to_two_ulp() {
        for (y, expected) in REFERENCE {
            let got = [f_7_6(y), f_4_3(y), f_11_6(y), f_2_3(y)];
            for (g, e) in got.iter().zip(expected) {
                assert!(
                    (g - e).abs() <= 2.0 * f64::EPSILON * e,
                    "y = {y}: {g} versus {e}"
                );
            }
        }
    }

    #[test]
    fn series_are_one_at_zero() {
        assert_eq!(f_7_6(0.0), 1.0);
        assert_eq!(f_4_3(0.0), 1.0);
        assert_eq!(f_11_6(0.0), 1.0);
        assert_eq!(f_2_3(0.0), 1.0);
    }
}
