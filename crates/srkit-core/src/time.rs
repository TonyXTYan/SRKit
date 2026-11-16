//! Exact timestamps and signed differences: whole ticks in a fixed-width unsigned integer.
//!
//! A timestamp counts ticks and does not know how long a tick is. The width is the type's limb count
//! ([`Time384`], [`Time512`], [`Time768`], [`Time1024`]); the tick's length belongs to a
//! [`crate::causal::Resolution`], and [`crate::causal::Clock`] holds a width and a resolution together. Methods
//! with `millis` in their name, and those without `_at`, read a tick as a millisecond, the default; the `_at`
//! methods take a resolution.
//!
//! Epochs are never converted to floating point. Subtract epochs first, then use
//! [`Timestamp::local_seconds_since`] only within its documented local range. Decimal text is the lossless
//! interchange format.

use bnum::BUint;
use core::{fmt, str::FromStr};

use crate::{causal::Resolution, Error};

/// Largest difference in ticks that the local float API accepts: beyond 2^53 a float no longer holds every
/// tick. That is 285,000 years of milliseconds and 104 days of nanoseconds.
pub const MAX_LOCAL_MILLIS: u64 = 1_u64 << 53;

/// Unsigned, fixed-width ticks since an application-defined epoch; a tick is a millisecond unless a
/// [`Resolution`] says otherwise.
///
/// `N` is the number of 64-bit limbs and must be positive: [`Time384`], [`Time512`], [`Time768`] and
/// [`Time1024`] are the widths that are tested and benchmarked.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Timestamp<const N: usize> {
    millis: BUint<N>,
}

/// The U384 millisecond timestamp prototype.
pub type Time384 = Timestamp<6>;
/// The U512 millisecond timestamp prototype.
pub type Time512 = Timestamp<8>;
/// Signed duration with a U384 magnitude.
pub type Duration384 = Duration<6>;
/// Signed duration with a U512 magnitude.
pub type Duration512 = Duration<8>;
/// A 768-bit timestamp. At one nanosecond a tick it reaches 4.9e214 years; at a picosecond, 4.9e211.
pub type Time768 = Timestamp<12>;
/// Signed duration with a U768 magnitude.
pub type Duration768 = Duration<12>;
/// A 1,024-bit timestamp (R71.4). At one nanosecond a tick it reaches 5.7e291 years; 512 bits at a nanosecond
/// reach 4.2e137 years, and 384 bits fall just short of 1e100 years.
pub type Time1024 = Timestamp<16>;
/// Signed duration with a U1024 magnitude.
pub type Duration1024 = Duration<16>;

/// A signed difference in ticks, with a full-width unsigned magnitude.
///
/// Negative zero is normalized to zero. Using sign and magnitude preserves
/// the difference between any two timestamps of the same width.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Duration<const N: usize> {
    negative: bool,
    magnitude: BUint<N>,
}

fn parse_unsigned<const N: usize>(text: &str) -> Result<BUint<N>, Error> {
    if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(Error::InvalidInput);
    }
    BUint::<N>::from_str_radix(text, 10).map_err(|_| Error::Overflow)
}

impl<const N: usize> Timestamp<N> {
    /// Epoch zero.
    pub const ZERO: Self = Self {
        millis: BUint::ZERO,
    };
    /// Maximum representable epoch; additions remain checked.
    pub const MAX: Self = Self { millis: BUint::MAX };

    /// Construct from a 64-bit count of ticks. The type counts ticks and does not know their length: the
    /// default tick is a millisecond, which the `*_millis` and `local_seconds` methods assume, and
    /// [`crate::causal::Resolution`] names other lengths for the exact causal checks.
    pub const fn from_ticks(ticks: u64) -> Self {
        Self::from_millis(ticks)
    }

    /// Construct an exact timestamp from a 64-bit count of milliseconds.
    pub const fn from_millis(millis: u64) -> Self {
        Self {
            millis: BUint::from_digit(millis),
        }
    }

    /// Construct from a full-width millisecond count.
    pub const fn from_raw_millis(millis: BUint<N>) -> Self {
        Self { millis }
    }

    /// Access the exact full-width count without conversion to floating point.
    pub const fn raw_millis(self) -> BUint<N> {
        self.millis
    }

    /// Construct from a full-width count of ticks, whatever their length.
    pub const fn from_raw_ticks(ticks: BUint<N>) -> Self {
        Self { millis: ticks }
    }

    /// The exact full-width count of ticks.
    pub const fn raw_ticks(self) -> BUint<N> {
        self.millis
    }

    /// Parse nonempty unsigned decimal digits, rejecting signs and whitespace.
    pub fn from_decimal(text: &str) -> Result<Self, Error> {
        Ok(Self {
            millis: parse_unsigned(text)?,
        })
    }

    /// Add a signed duration, rejecting wraparound and times before epoch zero.
    pub fn checked_add(self, duration: Duration<N>) -> Result<Self, Error> {
        let millis = if duration.negative {
            self.millis.checked_sub(duration.magnitude)
        } else {
            self.millis.checked_add(duration.magnitude)
        }
        .ok_or(Error::Overflow)?;
        Ok(Self { millis })
    }

    /// Subtract a signed duration without wrapping.
    pub fn checked_sub(self, duration: Duration<N>) -> Result<Self, Error> {
        self.checked_add(duration.negated())
    }

    /// Exact signed difference, including differences spanning the full width.
    pub fn signed_duration_since(self, origin: Self) -> Duration<N> {
        if self >= origin {
            Duration::from_parts(false, self.millis - origin.millis)
        } else {
            Duration::from_parts(true, origin.millis - self.millis)
        }
    }

    /// Exact nonnegative difference; an earlier timestamp is out of range.
    pub fn duration_since(self, origin: Self) -> Result<Duration<N>, Error> {
        if self < origin {
            return Err(Error::OutOfRange);
        }
        Ok(self.signed_duration_since(origin))
    }

    /// Subtract in integer arithmetic, then convert the local difference to s.
    ///
    /// Differences beyond ±2^53 milliseconds return `PrecisionLoss`. The
    /// integer millisecond value is exactly representable before division by
    /// 1000; seconds receive the usual binary64 rounding. No large-epoch float
    /// subtraction occurs.
    pub fn local_seconds_since(self, origin: Self) -> Result<f64, Error> {
        self.signed_duration_since(origin).local_seconds()
    }

    /// Exact integer phase for a positive millisecond period, even at huge epochs.
    pub fn modulo_millis(self, period_millis: u64) -> Result<u64, Error> {
        if period_millis == 0 {
            return Err(Error::InvalidInput);
        }
        let phase = self.millis % BUint::<N>::from_digit(period_millis);
        Ok(phase.digits()[0])
    }

    /// Exact integer phase in a positive period of whole ticks, even at huge epochs.
    pub fn modulo_ticks(self, period_ticks: u64) -> Result<u64, Error> {
        self.modulo_millis(period_ticks)
    }

    /// Exact integer phase in a period of any length. A period of 64 bits is 584 million years of milliseconds
    /// but seven months of picoseconds, so a fine tick needs this form for a period as long as a year.
    pub fn modulo(self, period_ticks: BUint<N>) -> Result<BUint<N>, Error> {
        if period_ticks == BUint::ZERO {
            return Err(Error::InvalidInput);
        }
        Ok(self.millis % period_ticks)
    }

    /// As [`Timestamp::local_seconds_since`], for ticks of the length `res` gives. Differences beyond ±2^53
    /// ticks return `PrecisionLoss`.
    pub fn local_seconds_since_at(self, res: Resolution, origin: Self) -> Result<f64, Error> {
        self.signed_duration_since(origin).local_seconds_at(res)
    }

    /// The exact phase in a period of whole ticks, in seconds, for ticks of the length `res` gives.
    pub fn phase_seconds_at(self, res: Resolution, period_ticks: u64) -> Result<f64, Error> {
        let phase = self.modulo_ticks(period_ticks)?;
        if phase > MAX_LOCAL_MILLIS {
            return Err(Error::PrecisionLoss);
        }
        Ok(phase as f64 / res.ticks_per_second() as f64)
    }

    /// Convert the exact modular phase to local seconds, after reduction.
    pub fn phase_seconds(self, period_millis: u64) -> Result<f64, Error> {
        let phase = self.modulo_millis(period_millis)?;
        if phase > MAX_LOCAL_MILLIS {
            return Err(Error::PrecisionLoss);
        }
        Ok(phase as f64 / 1000.0)
    }

    /// First causally allowed tick for a signal crossing a fixed integer-metre
    /// distance from this emission timestamp.
    ///
    /// This is an exact stationary-distance predicate, not a claim that a
    /// floating-point worldline intersection has exact geometry. The SI speed
    /// of light is exact; arrival is rounded upward to the next millisecond.
    pub fn stationary_arrival_metres(self, distance_metres: BUint<N>) -> Result<Self, Error> {
        self.checked_add(Duration::light_travel_time_metres(distance_metres)?)
    }

    /// Whether this tick is causally late enough for the given fixed distance.
    ///
    /// Uses integer arithmetic throughout and returns false before emission.
    /// Unlike constructing an arrival timestamp, this can answer false even
    /// when emission plus the travel time would exceed the timestamp range.
    pub fn stationary_signal_arrived(
        self,
        emission: Self,
        distance_metres: BUint<N>,
    ) -> Result<bool, Error> {
        if self < emission {
            return Ok(false);
        }
        let delay = Duration::<N>::light_travel_time_metres(distance_metres)?;
        Ok(self.signed_duration_since(emission).magnitude >= delay.magnitude)
    }
}

impl<const N: usize> Duration<N> {
    /// Zero has no negative sign.
    pub const ZERO: Self = Self {
        negative: false,
        magnitude: BUint::ZERO,
    };

    /// Exact causal travel time across a fixed nonnegative integer-metre
    /// distance, rounded up to integer milliseconds.
    ///
    /// Computes `ceil(distance * 1000 / 299792458)` without overflowing an
    /// intermediate `distance * 1000`, even at full-width distances. If
    /// `distance = q*c + r`, the result is `1000*q + ceil(1000*r/c)`.
    pub fn light_travel_time_metres(distance_metres: BUint<N>) -> Result<Self, Error> {
        const C_METRES_PER_SECOND: u64 = 299_792_458;
        let c = BUint::<N>::from_digit(C_METRES_PER_SECOND);
        let whole_seconds = distance_metres / c;
        let remaining_metres = (distance_metres % c).digits()[0];
        let residual_millis = (remaining_metres * 1000).div_ceil(C_METRES_PER_SECOND);
        let magnitude = whole_seconds
            .checked_mul(BUint::from_digit(1000))
            .and_then(|millis| millis.checked_add(BUint::from_digit(residual_millis)))
            .ok_or(Error::Overflow)?;
        Ok(Self::from_parts(false, magnitude))
    }

    /// Construct a signed duration from a 64-bit millisecond count.
    pub fn from_millis(millis: i64) -> Self {
        Self::from_parts(millis < 0, BUint::from_digit(millis.unsigned_abs()))
    }

    /// Construct a nonnegative duration from a 64-bit count of ticks (see [`Timestamp::from_ticks`]).
    pub const fn from_unsigned_ticks(ticks: u64) -> Self {
        Self::from_unsigned_millis(ticks)
    }

    /// Construct a nonnegative duration from a full-width count of ticks.
    pub const fn from_raw_ticks(ticks: BUint<N>) -> Self {
        Self {
            negative: false,
            magnitude: ticks,
        }
    }

    /// Construct a nonnegative duration from a 64-bit millisecond count.
    pub const fn from_unsigned_millis(millis: u64) -> Self {
        Self {
            negative: false,
            magnitude: BUint::from_digit(millis),
        }
    }

    /// Construct a full-width duration, normalizing negative zero.
    pub fn from_parts(negative: bool, magnitude: BUint<N>) -> Self {
        Self {
            negative: negative && magnitude != BUint::ZERO,
            magnitude,
        }
    }

    /// Parse unsigned decimal digits with an optional leading minus sign.
    pub fn from_decimal(text: &str) -> Result<Self, Error> {
        let (negative, digits) = match text.strip_prefix('-') {
            Some(digits) => (true, digits),
            None => (false, text),
        };
        Ok(Self::from_parts(negative, parse_unsigned(digits)?))
    }

    /// Whether this is strictly less than zero.
    pub const fn is_negative(self) -> bool {
        self.negative
    }

    /// The exact unsigned magnitude in milliseconds.
    pub const fn magnitude_millis(self) -> BUint<N> {
        self.magnitude
    }

    /// The exact unsigned magnitude in ticks, whatever their length.
    pub const fn magnitude_ticks(self) -> BUint<N> {
        self.magnitude
    }

    /// Construct a signed duration from a 64-bit count of ticks.
    pub fn from_ticks(ticks: i64) -> Self {
        Self::from_millis(ticks)
    }

    /// Local seconds for ticks of the length `res` gives. Magnitudes beyond 2^53 ticks return `PrecisionLoss`;
    /// below that the count is exact before the one division.
    pub fn local_seconds_at(self, res: Resolution) -> Result<f64, Error> {
        if self.magnitude > BUint::<N>::from_digit(MAX_LOCAL_MILLIS) {
            return Err(Error::PrecisionLoss);
        }
        let seconds = self.magnitude.digits()[0] as f64 / res.ticks_per_second() as f64;
        Ok(if self.negative { -seconds } else { seconds })
    }

    /// Negation cannot overflow a sign-and-magnitude duration.
    pub fn negated(self) -> Self {
        Self::from_parts(!self.negative, self.magnitude)
    }

    /// Signed addition, failing only when the resulting magnitude overflows.
    pub fn checked_add(self, other: Self) -> Result<Self, Error> {
        if self.negative == other.negative {
            let magnitude = self
                .magnitude
                .checked_add(other.magnitude)
                .ok_or(Error::Overflow)?;
            Ok(Self::from_parts(self.negative, magnitude))
        } else if self.magnitude >= other.magnitude {
            Ok(Self::from_parts(
                self.negative,
                self.magnitude - other.magnitude,
            ))
        } else {
            Ok(Self::from_parts(
                other.negative,
                other.magnitude - self.magnitude,
            ))
        }
    }

    /// Signed subtraction without wrapping.
    pub fn checked_sub(self, other: Self) -> Result<Self, Error> {
        self.checked_add(other.negated())
    }

    /// Local seconds; see [`Timestamp::local_seconds_since`] for precision limits.
    pub fn local_seconds(self) -> Result<f64, Error> {
        if self.magnitude > BUint::<N>::from_digit(MAX_LOCAL_MILLIS) {
            return Err(Error::PrecisionLoss);
        }
        let seconds = self.magnitude.digits()[0] as f64 / 1000.0;
        Ok(if self.negative { -seconds } else { seconds })
    }
}

impl<const N: usize> fmt::Display for Timestamp<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.millis.fmt(f)
    }
}

impl<const N: usize> fmt::Display for Duration<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.negative {
            f.write_str("-")?;
        }
        self.magnitude.fmt(f)
    }
}

impl<const N: usize> FromStr for Timestamp<N> {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::from_decimal(text)
    }
}

impl<const N: usize> FromStr for Duration<N> {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::from_decimal(text)
    }
}
