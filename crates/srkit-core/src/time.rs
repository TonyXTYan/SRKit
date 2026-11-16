//! Exact millisecond timestamps and signed differences.
//!
//! Epochs are never converted to floating point. Subtract epochs first, then
//! use [`Timestamp::local_seconds_since`] only within its documented local
//! range. Both widths remain available while the project's width comparison
//! is open. Decimal text is the lossless interchange format.

use bnum::BUint;
use core::{fmt, str::FromStr};

use crate::Error;

/// Largest integer millisecond difference accepted by the local float API.
pub const MAX_LOCAL_MILLIS: u64 = 1_u64 << 53;

/// Unsigned, fixed-width milliseconds since an application-defined epoch.
///
/// `N` is the number of 64-bit limbs and must be positive. Use [`Time384`]
/// or [`Time512`] unless deliberately testing another width.
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

/// A signed millisecond difference, with a full-width unsigned magnitude.
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
