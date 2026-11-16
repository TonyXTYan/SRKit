//! Axis ticks at "nice" steps, with labels built from integers so that no tick reads `0.30000000000000004`.

/// Powers of ten that binary64 holds exactly.
const POW10: [f64; 23] = [
    1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16,
    1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
];

fn pow10(exponent: i32) -> f64 {
    match usize::try_from(exponent) {
        Ok(i) if i < POW10.len() => POW10[i],
        _ => libm::pow(10.0, f64::from(exponent)),
    }
}

/// `value × 10^exponent`, dividing for a negative exponent so that `3 × 10^-1` is the double nearest 0.3.
/// Below 10^-300 the division is done in two parts, since the divisor alone would overflow.
fn scale10(value: f64, exponent: i32) -> f64 {
    if exponent >= 0 {
        value * pow10(exponent)
    } else if exponent >= -300 {
        value / pow10(-exponent)
    } else {
        value / 1e300 / pow10(-exponent - 300)
    }
}

/// A tick step `mantissa × 10^exponent` with a mantissa of 1, 2 or 5.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Step {
    pub mantissa: i64,
    pub exponent: i32,
}

impl Step {
    pub fn value(self) -> f64 {
        scale10(self.mantissa as f64, self.exponent)
    }
    /// The next larger nice step.
    pub fn next(self) -> Self {
        match self.mantissa {
            1 => Self {
                mantissa: 2,
                ..self
            },
            2 => Self {
                mantissa: 5,
                ..self
            },
            _ => Self {
                mantissa: 1,
                exponent: self.exponent + 1,
            },
        }
    }
    /// The smallest nice step that is at least `raw`, which must be positive and finite.
    pub fn at_least(raw: f64) -> Self {
        // The logarithm may be one off at a power of ten; start one decade low and walk up. The clamp keeps
        // the exponent inside the range of binary64 whatever comes in.
        let decade = libm::floor(libm::log10(raw));
        let exponent = if decade.is_finite() {
            decade.clamp(-330.0, 330.0) as i32 - 1
        } else {
            0
        };
        let mut step = Self {
            mantissa: 1,
            exponent,
        };
        // Three steps to a decade: 32 covers the slack with room to spare, and bounds the loop for any input.
        for _ in 0..32 {
            if step.value() >= raw {
                break;
            }
            step = step.next();
        }
        step
    }
    /// The step as text: a decimal string for steps of ordinary size (`0.5`, `2`, `1000`), a power of
    /// ten otherwise (`5e-9`, `2e12`).
    pub fn text(self) -> String {
        if (-4..=5).contains(&self.exponent) {
            decimal(self.mantissa, self.exponent)
        } else {
            format!("{}e{}", self.mantissa, self.exponent)
        }
    }
}

/// `n × 10^exponent` as a decimal string with no rounding: digits are shifted, never computed in floating
/// point. Negative numbers take the minus sign U+2212.
pub(crate) fn decimal(n: i64, exponent: i32) -> String {
    let mut digits = n.unsigned_abs().to_string();
    if exponent >= 0 {
        if n != 0 {
            digits.extend(std::iter::repeat('0').take(exponent as usize));
        }
    } else {
        let places = exponent.unsigned_abs() as usize;
        if digits.len() <= places {
            digits = format!("{digits:0>width$}", width = places + 1);
        }
        digits.insert(digits.len() - places, '.');
    }
    if n < 0 {
        digits.insert(0, '\u{2212}');
    }
    digits
}

/// The power of a thousand to take out of numbers whose largest has the decimal exponent `top`: none for
/// 0.001 to 99999, which read well as plain decimals.
fn thousands(top: i32) -> i32 {
    if (-3..5).contains(&top) {
        0
    } else {
        3 * top.div_euclid(3)
    }
}

/// The readings of a clock at the multiples of a step, written exactly: the step is taken as an integer of at
/// most six digits times a power of ten, and each reading is that integer multiplied out, so 3 × 0.1 reads
/// `0.3`. Every reading has the same number of decimals. Where the largest is below 0.001 or at least
/// 100000, all are written with one shared power of ten: `0.5e-6`, `1.0e-6`, `1.5e-6`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Readings {
    mantissa: i64,
    exponent: i32,
    /// The power of ten written after each reading. Zero for plain decimals.
    shared: i32,
}

impl Readings {
    /// Readings at the multiples `first..=last` of `step`, which must be positive and finite.
    pub fn new(step: f64, first: i64, last: i64) -> Self {
        // Six significant digits of the step as an integer, then without its trailing zeros.
        let decade = libm::floor(libm::log10(step));
        let mut exponent = if decade.is_finite() {
            decade.clamp(-330.0, 330.0) as i32 - 5
        } else {
            0
        };
        // `as` saturates; a step that is not a number gives a mantissa of zero and readings of zero.
        let mut mantissa = scale10(step, -exponent).round() as i64;
        while mantissa != 0 && mantissa % 10 == 0 {
            mantissa /= 10;
            exponent += 1;
        }
        let largest = first
            .unsigned_abs()
            .max(last.unsigned_abs())
            .saturating_mul(mantissa.unsigned_abs());
        let top = if largest == 0 {
            0
        } else {
            largest.ilog10() as i32 + exponent
        };
        Self {
            mantissa,
            exponent,
            shared: thousands(top),
        }
    }
    /// The reading at multiple `k` of the step.
    pub fn text(&self, k: i64) -> String {
        let digits = decimal(k.saturating_mul(self.mantissa), self.exponent - self.shared);
        if self.shared == 0 {
            digits
        } else {
            format!("{digits}e{}", self.shared)
        }
    }
}

/// The ticks of one axis: every multiple of a nice step inside a range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Ticks {
    pub step: Step,
    /// First and last multiples of the step. Empty when `first > last`.
    pub first: i64,
    pub last: i64,
    /// Power of ten taken out of the labels and shown once beside the axis title. Zero for plain labels.
    pub exponent: i32,
}

impl Ticks {
    /// The multiples of `step` in `[min, max]`. A bound that is a multiple to within rounding counts as inside.
    pub fn new(min: f64, max: f64, step: Step) -> Self {
        let size = step.value();
        let (lo, hi) = (min / size, max / size);
        let slack = |q: f64| (q.abs() * 1e-12).max(1e-9);
        // `as` saturates, so a quotient beyond the integers cannot wrap.
        let first = (lo - slack(lo)).ceil() as i64;
        let last = (hi + slack(hi)).floor() as i64;
        let largest = first
            .unsigned_abs()
            .max(last.unsigned_abs())
            .saturating_mul(step.mantissa.unsigned_abs());
        // Decimal exponent of the largest label: digits - 1 + step exponent.
        let top = if first > last || largest == 0 {
            0
        } else {
            largest.ilog10() as i32 + step.exponent
        };
        Self {
            step,
            first,
            last,
            exponent: thousands(top),
        }
    }
    pub fn count(&self) -> usize {
        if self.first > self.last {
            0
        } else {
            usize::try_from(self.last.saturating_sub(self.first).saturating_add(1))
                .unwrap_or(usize::MAX)
        }
    }
    /// Tick indices, first to last.
    pub fn indices(&self) -> impl Iterator<Item = i64> {
        self.first..=self.last
    }
    /// The coordinate of tick `k`.
    pub fn value(&self, k: i64) -> f64 {
        scale10(
            k.saturating_mul(self.step.mantissa) as f64,
            self.step.exponent,
        )
    }
    /// The label of tick `k`, in units of `10^exponent`.
    pub fn label(&self, k: i64) -> String {
        decimal(
            k.saturating_mul(self.step.mantissa),
            self.step.exponent - self.exponent,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steps_are_one_two_five() {
        let cases = [
            (0.7, 1, 0),
            (1.0, 1, 0),
            (1.0001, 2, 0),
            (2.0, 2, 0),
            (3.0, 5, 0),
            (5.0, 5, 0),
            (6.0, 1, 1),
            (0.03, 5, -2),
            (0.1, 1, -1),
            (1234.0, 2, 3),
            (1e-7, 1, -7),
            (4e22, 5, 22),
        ];
        for (raw, mantissa, exponent) in cases {
            let step = Step::at_least(raw);
            assert_eq!(
                (step.mantissa, step.exponent),
                (mantissa, exponent),
                "{raw}"
            );
            assert!(step.value() >= raw);
        }
        let mut step = Step::at_least(0.01);
        let mut seen = Vec::new();
        for _ in 0..7 {
            seen.push(step.text());
            step = step.next();
        }
        assert_eq!(seen, ["0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1"]);
        assert_eq!(Step::at_least(4e-9).text(), "5e-9");
        assert_eq!(Step::at_least(1.5e12).text(), "2e12");
        assert_eq!(Step::at_least(100_000.0).text(), "100000");
        // numbers at the ends of the range of binary64 still give a step
        for raw in [5e-324, f64::MIN_POSITIVE, 1e300, f64::MAX] {
            let step = Step::at_least(raw);
            assert!(
                step.value() > 0.0 && (step.value() >= raw || raw > 1e308),
                "{raw:e}"
            );
        }
    }

    #[test]
    fn decimal_strings_are_exact() {
        assert_eq!(decimal(3, -1), "0.3");
        assert_eq!(decimal(-3, -1), "\u{2212}0.3");
        assert_eq!(decimal(15, -1), "1.5");
        assert_eq!(decimal(0, -2), "0.00");
        assert_eq!(decimal(0, 3), "0");
        assert_eq!(decimal(25, 2), "2500");
        assert_eq!(decimal(-5, -3), "\u{2212}0.005");
        assert_eq!(decimal(120, -2), "1.20");
    }

    fn labels(min: f64, max: f64, raw: f64) -> (Vec<String>, i32) {
        let ticks = Ticks::new(min, max, Step::at_least(raw));
        (
            ticks.indices().map(|k| ticks.label(k)).collect(),
            ticks.exponent,
        )
    }

    #[test]
    fn ticks_cover_the_range_with_exact_labels() {
        assert_eq!(
            labels(-0.3, 0.3, 0.1),
            (
                [
                    "\u{2212}0.3",
                    "\u{2212}0.2",
                    "\u{2212}0.1",
                    "0.0",
                    "0.1",
                    "0.2",
                    "0.3"
                ]
                .map(String::from)
                .to_vec(),
                0
            )
        );
        assert_eq!(
            labels(-0.5, 12.3, 2.0),
            (
                ["0", "2", "4", "6", "8", "10", "12"]
                    .map(String::from)
                    .to_vec(),
                0
            )
        );
        // the tick values are the doubles nearest the decimals
        let ticks = Ticks::new(0.0, 1.0, Step::at_least(0.1));
        assert_eq!(ticks.value(3), 0.3);
        assert_eq!(ticks.value(7), 0.7);
        assert_eq!(ticks.count(), 11);
    }

    #[test]
    fn large_and_small_ranges_take_out_a_power_of_a_thousand() {
        assert_eq!(
            labels(0.0, 1.4e8, 5e7),
            (["0", "50", "100"].map(String::from).to_vec(), 6)
        );
        assert_eq!(
            labels(0.0, 3e-5, 1e-5),
            (["0", "10", "20", "30"].map(String::from).to_vec(), -6)
        );
        assert_eq!(
            labels(0.0, 250_000.0, 100_000.0),
            (["0", "100", "200"].map(String::from).to_vec(), 3)
        );
        assert_eq!(
            labels(0.0, 0.004, 0.002),
            (["0.000", "0.002", "0.004"].map(String::from).to_vec(), 0)
        );
        assert_eq!(labels(0.0, 99_999.0, 50_000.0).1, 0);
    }

    #[test]
    fn clock_readings_are_exact_at_any_scale() {
        let texts = |step: f64, first: i64, last: i64| -> Vec<String> {
            let readings = Readings::new(step, first, last);
            (first..=last).map(|k| readings.text(k)).collect()
        };
        assert_eq!(texts(1.0, 0, 3), ["0", "1", "2", "3"]);
        assert_eq!(texts(0.5, 0, 3), ["0.0", "0.5", "1.0", "1.5"]);
        assert_eq!(
            texts(0.75, -2, 2),
            ["\u{2212}1.50", "\u{2212}0.75", "0.00", "0.75", "1.50"]
        );
        // 0.1 is not a binary fraction, but its multiples are written as decimal ones
        assert_eq!(texts(0.1, 2, 4), ["0.2", "0.3", "0.4"]);
        assert_eq!(texts(20.0, 0, 2), ["0", "20", "40"]);
        assert_eq!(texts(0.002, 1, 3), ["0.002", "0.004", "0.006"]);
        // small and large readings share one power of a thousand
        assert_eq!(
            texts(5e-7, 0, 4),
            ["0.0e-6", "0.5e-6", "1.0e-6", "1.5e-6", "2.0e-6"]
        );
        assert_eq!(texts(5e6, 1, 3), ["5e6", "10e6", "15e6"]);
        assert_eq!(texts(2.5e-10, 3, 4), ["0.75e-9", "1.00e-9"]);
        // a step with more than six digits is rounded to six
        assert_eq!(texts(1.0 / 3.0, 1, 3), ["0.333333", "0.666666", "0.999999"]);
        // steps at the ends of the range of binary64 give text, not a panic
        for step in [
            5e-324,
            f64::MIN_POSITIVE,
            1e300,
            f64::MAX,
            f64::NAN,
            f64::INFINITY,
            0.0,
        ] {
            assert!(!Readings::new(step, i64::MIN, i64::MAX)
                .text(i64::MAX)
                .is_empty());
        }
    }

    #[test]
    fn an_empty_range_has_no_ticks() {
        let ticks = Ticks::new(0.1, 0.9, Step::at_least(1.0));
        assert_eq!(ticks.count(), 0);
        assert_eq!(ticks.indices().count(), 0);
        assert_eq!(ticks.exponent, 0);
    }
}
