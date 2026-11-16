//! Exact causal decisions on stored integer times.
//!
//! Geometry in this crate is floating point, which is right for drawing and wrong for deciding. "Has this
//! signal arrived by tick n?" and "is event B inside A's light cone?" are answered here in integer arithmetic:
//! times are whole milliseconds ([`crate::time::Timestamp`]) and positions are whole metres on a grid
//! ([`GridVec`]). The speed of light is exactly 299,792,458 m/s, so
//!
//! ```text
//! arrived  <=>  c * dt_ms >= 1000 * |dx_m|  <=>  (c * dt_ms)² >= 10⁶ * (dx² + dy² + dz²)
//! ```
//!
//! holds or fails with no rounding at all. This is the "exact predicates, inexact constructions" split of
//! computational geometry: compute where things are in floating point, round the endpoints to the grid, then
//! decide exactly.
//!
//! The grid is one metre and the tick one millisecond, in which light covers about 300 km: rounding a position
//! to the grid moves an arrival by at most one part in 300,000 of a tick.

use bnum::BUint;

use crate::{
    time::{Duration, Timestamp},
    Error,
};

/// Speed of light in metres per second, exactly.
pub const C_METRES_PER_SECOND: u64 = 299_792_458;

/// The working width for squared quantities: 512 bits holds `10⁶ * 3 * (2^127)²` with room to spare.
type Wide = BUint<8>;

/// A displacement in whole metres. `i128` metres span 1.8 × 10²² light-years each way.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct GridVec {
    pub x: i128,
    pub y: i128,
    pub z: i128,
}

impl GridVec {
    pub const ZERO: Self = Self::new(0, 0, 0);
    pub const fn new(x: i128, y: i128, z: i128) -> Self {
        Self { x, y, z }
    }
    /// Round a floating-point displacement in metres to the grid (to nearest, ties away from zero).
    pub fn from_metres(x: f64, y: f64, z: f64) -> Result<Self, Error> {
        let conv = |v: f64| {
            if !v.is_finite() || v.abs() >= 1.7e38 {
                Err(Error::Overflow)
            } else {
                Ok(libm::round(v) as i128)
            }
        };
        Ok(Self::new(conv(x)?, conv(y)?, conv(z)?))
    }
    /// The displacement from `self` to `other`.
    pub fn to(self, other: Self) -> Result<Self, Error> {
        Ok(Self::new(
            other.x.checked_sub(self.x).ok_or(Error::Overflow)?,
            other.y.checked_sub(self.y).ok_or(Error::Overflow)?,
            other.z.checked_sub(self.z).ok_or(Error::Overflow)?,
        ))
    }
    /// Squared length in square metres, exactly.
    fn norm_squared(self) -> Wide {
        let sq = |v: i128| {
            let w = wide_from_u128(v.unsigned_abs());
            w * w
        };
        sq(self.x) + sq(self.y) + sq(self.z)
    }
    /// The sum of the absolute components, which is never less than the length.
    fn l1(self) -> Wide {
        wide_from_u128(self.x.unsigned_abs())
            + wide_from_u128(self.y.unsigned_abs())
            + wide_from_u128(self.z.unsigned_abs())
    }
}

fn wide_from_u128(v: u128) -> Wide {
    let mut digits = [0u64; 8];
    digits[0] = v as u64;
    digits[1] = (v >> 64) as u64;
    Wide::from_digits(digits)
}

/// Copy a magnitude into the working width, or `None` if it does not fit in 512 bits.
fn widen<const N: usize>(v: BUint<N>) -> Option<Wide> {
    let src = v.digits();
    let mut digits = [0u64; 8];
    for (i, d) in src.iter().enumerate() {
        if i < 8 {
            digits[i] = *d;
        } else if *d != 0 {
            return None;
        }
    }
    Some(Wide::from_digits(digits))
}

fn narrow<const N: usize>(v: Wide) -> Option<BUint<N>> {
    let src = v.digits();
    let mut digits = [0u64; N];
    for (i, d) in src.iter().enumerate() {
        if i < N {
            digits[i] = *d;
        } else if *d != 0 {
            return None;
        }
    }
    Some(BUint::from_digits(digits))
}

/// The smallest integer whose square is at least `v`.
fn isqrt_ceil(v: Wide) -> Wide {
    if v == Wide::ZERO {
        return Wide::ZERO;
    }
    // Newton's iteration from above converges to floor(sqrt(v)).
    let bits = v.bits();
    let mut x = Wide::ONE << bits.div_ceil(2);
    loop {
        let y = (x + v / x) >> 1;
        if y >= x {
            break;
        }
        x = y;
    }
    if x * x == v {
        x
    } else {
        x + Wide::ONE
    }
}

/// How two events on the grid are related, exactly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Causal {
    /// The later event is strictly inside the earlier one's light cone.
    Timelike,
    /// Light from the earlier event reaches the later one exactly.
    Null,
    /// Neither event can influence the other.
    Spacelike,
}

/// What one tick of the clock and one step of the grid are.
///
/// The game stores milliseconds on a metre grid ([`Resolution::MILLISECOND_METRE`]), and the functions without
/// `_at` in their name assume it. The `_at` functions take any resolution, for trying finer clocks (R71.4).
///
/// The grid must be much finer than the distance light covers in a tick, or rounding a position to the grid
/// moves an arrival by whole ticks: see [`Resolution::light_steps_per_tick`]. A nanosecond tick on a metre grid
/// fails that (light covers 0.3 m in a nanosecond); a nanosecond needs a millimetre grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Resolution {
    ticks_per_second: u64,
    grid_per_metre: u64,
}

impl Resolution {
    /// Milliseconds and metres: light crosses 299,792 grid steps in a tick. `i128` metres span 1.8e22 ly.
    pub const MILLISECOND_METRE: Self = Self {
        ticks_per_second: 1_000,
        grid_per_metre: 1,
    };
    /// Microseconds and metres: light crosses 300 grid steps in a tick.
    pub const MICROSECOND_METRE: Self = Self {
        ticks_per_second: 1_000_000,
        grid_per_metre: 1,
    };
    /// Nanoseconds and millimetres: light crosses 300 grid steps in a tick. `i128` millimetres span 1.8e19 ly.
    pub const NANOSECOND_MILLIMETRE: Self = Self {
        ticks_per_second: 1_000_000_000,
        grid_per_metre: 1_000,
    };

    /// A resolution of `ticks_per_second` ticks to the second and `grid_per_metre` grid steps to the metre.
    /// Both must be at least 1; the tick no finer than a picosecond and the grid no finer than a micrometre,
    /// which keeps every intermediate product within 512 bits.
    pub const fn new(ticks_per_second: u64, grid_per_metre: u64) -> Result<Self, Error> {
        if ticks_per_second == 0
            || grid_per_metre == 0
            || ticks_per_second > 1_000_000_000_000
            || grid_per_metre > 1_000_000
        {
            return Err(Error::InvalidInput);
        }
        Ok(Self {
            ticks_per_second,
            grid_per_metre,
        })
    }
    pub const fn ticks_per_second(self) -> u64 {
        self.ticks_per_second
    }
    pub const fn grid_per_metre(self) -> u64 {
        self.grid_per_metre
    }
    /// How many grid steps light covers in one tick. Rounding a position to the grid moves an arrival by up to
    /// one part in this many of a tick, so it should be large (hundreds or more).
    pub fn light_steps_per_tick(self) -> f64 {
        C_METRES_PER_SECOND as f64 * self.grid_per_metre as f64 / self.ticks_per_second as f64
    }
}

/// Compare `c * |dt|` with `|dx|` exactly, in ticks and grid steps: the sign of
/// `(c g dt)² - (f |dx|)²`, where `f` is ticks per second and `g` grid steps per metre.
fn compare<const N: usize>(res: Resolution, dt: BUint<N>, dx: GridVec) -> core::cmp::Ordering {
    let speed = Wide::from_digit(C_METRES_PER_SECOND) * Wide::from_digit(res.grid_per_metre);
    let rate = Wide::from_digit(res.ticks_per_second);
    // A wait longer than light needs to cross the L1 distance is certainly long enough, and anything shorter
    // fits the working width with room for its square.
    let bound = dx.l1() * rate;
    let lhs = match widen(dt).and_then(|w| w.checked_mul(speed)) {
        Some(lhs) if lhs <= bound => lhs,
        _ => return core::cmp::Ordering::Greater,
    };
    let rhs = dx.norm_squared() * rate * rate;
    (lhs * lhs).cmp(&rhs)
}

/// Classify the separation between two events `dt` apart in time and `dx` apart in space, at a resolution.
pub fn classify_at<const N: usize>(res: Resolution, dt: Duration<N>, dx: GridVec) -> Causal {
    match compare(res, dt.magnitude_millis(), dx) {
        core::cmp::Ordering::Greater => Causal::Timelike,
        core::cmp::Ordering::Equal => Causal::Null,
        core::cmp::Ordering::Less => Causal::Spacelike,
    }
}

/// Classify the separation between two events `dt` milliseconds apart in time and `dx` metres apart in space.
pub fn classify<const N: usize>(dt: Duration<N>, dx: GridVec) -> Causal {
    classify_at(Resolution::MILLISECOND_METRE, dt, dx)
}

/// Whether a signal emitted at `emission` has reached a point `dx` away by the tick `now`, at a resolution.
/// Exact.
pub fn arrived_at<const N: usize>(
    res: Resolution,
    emission: Timestamp<N>,
    now: Timestamp<N>,
    dx: GridVec,
) -> bool {
    if now < emission {
        return false;
    }
    let dt = now.signed_duration_since(emission).magnitude_millis();
    compare(res, dt, dx) != core::cmp::Ordering::Less
}

/// Whether a signal emitted at `emission` has reached a point `dx` metres away by the millisecond `now`. Exact.
pub fn arrived<const N: usize>(emission: Timestamp<N>, now: Timestamp<N>, dx: GridVec) -> bool {
    arrived_at(Resolution::MILLISECOND_METRE, emission, now, dx)
}

/// The light-travel time over `dx`, rounded up to whole ticks: the smallest `n` with `c g n >= f |dx|`. Exact,
/// by an integer square root.
pub fn light_ticks_at<const N: usize>(res: Resolution, dx: GridVec) -> Result<Duration<N>, Error> {
    // c g n >= sqrt(f² |dx|²) and c g n is an integer, so c g n >= ceil(sqrt(...)).
    let rate = Wide::from_digit(res.ticks_per_second);
    let s = isqrt_ceil(dx.norm_squared() * rate * rate);
    let speed = Wide::from_digit(C_METRES_PER_SECOND) * Wide::from_digit(res.grid_per_metre);
    let n = s.div_ceil(speed);
    let magnitude = narrow::<N>(n).ok_or(Error::Overflow)?;
    Ok(Duration::from_parts(false, magnitude))
}

/// The light-travel time over `dx` metres, rounded up to whole milliseconds.
pub fn light_ticks<const N: usize>(dx: GridVec) -> Result<Duration<N>, Error> {
    light_ticks_at(Resolution::MILLISECOND_METRE, dx)
}

/// The first tick at which a signal emitted at `emission` has reached a fixed point `dx` away, at a resolution.
pub fn first_arrival_at<const N: usize>(
    res: Resolution,
    emission: Timestamp<N>,
    dx: GridVec,
) -> Result<Timestamp<N>, Error> {
    emission.checked_add(light_ticks_at::<N>(res, dx)?)
}

/// The first millisecond at which a signal emitted at `emission` has reached a fixed point `dx` metres away.
pub fn first_arrival<const N: usize>(
    emission: Timestamp<N>,
    dx: GridVec,
) -> Result<Timestamp<N>, Error> {
    first_arrival_at(Resolution::MILLISECOND_METRE, emission, dx)
}

/// The first tick at which a signal reaches a **moving** receiver, settled exactly from a floating-point
/// estimate.
///
/// `displacement_at(tick)` gives the receiver's position relative to the emission point at that tick, on the
/// grid (evaluate the receiver's worldline in floating point and round). Starting from `estimate`, the search
/// steps back while the signal has already arrived and forward until it has, so the result is the first tick
/// that passes the exact check `c dt >= |dx|`: it can never be a tick early. `guard` bounds the number of steps
/// each way; an estimate further off than that returns [`Error::PrecisionLoss`].
///
/// For a receiver slower than light the check, once true, stays true, so the answer is unique up to the
/// rounding of positions to the grid.
pub fn settle_arrival<const N: usize, F>(
    emission: Timestamp<N>,
    estimate: Timestamp<N>,
    guard: u32,
    displacement_at: F,
) -> Result<Timestamp<N>, Error>
where
    F: FnMut(Timestamp<N>) -> Result<GridVec, Error>,
{
    settle_arrival_at(
        Resolution::MILLISECOND_METRE,
        emission,
        estimate,
        guard,
        displacement_at,
    )
}

/// As [`settle_arrival`], at any resolution: ticks and grid steps are those of `res`.
pub fn settle_arrival_at<const N: usize, F>(
    res: Resolution,
    emission: Timestamp<N>,
    estimate: Timestamp<N>,
    guard: u32,
    mut displacement_at: F,
) -> Result<Timestamp<N>, Error>
where
    F: FnMut(Timestamp<N>) -> Result<GridVec, Error>,
{
    let one = Duration::<N>::from_unsigned_millis(1);
    let mut tick = if estimate < emission {
        emission
    } else {
        estimate
    };
    if arrived_at(res, emission, tick, displacement_at(tick)?) {
        for _ in 0..guard {
            if tick == emission {
                return Ok(tick);
            }
            let earlier = tick.checked_sub(one)?;
            if !arrived_at(res, emission, earlier, displacement_at(earlier)?) {
                return Ok(tick);
            }
            tick = earlier;
        }
    } else {
        for _ in 0..guard {
            tick = tick.checked_add(one)?;
            if arrived_at(res, emission, tick, displacement_at(tick)?) {
                return Ok(tick);
            }
        }
    }
    Err(Error::PrecisionLoss)
}
