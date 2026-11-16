//! `srkit.causal`: exact causal decisions on stored integer times.
//!
//! Times are whole ticks (`Time384`, `Time512`, `Time768`, `Time1024`) and displacements are whole grid steps
//! as Python `int`s, checked against the `i128` range. No float is involved in any decision.
//!
//! `Clock` holds a width, a tick and a grid spacing. The module-level functions are the same checks at the
//! default resolution, a millisecond tick on a metre grid, and take the width from the times they are given.
//! Both go through the functions below that are generic over the width; `at_width!` picks the width.

use bnum::BUint;
use pyo3::{
    basic::CompareOp,
    exceptions::{PyOverflowError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyInt, PyString, PyType},
};
use srkit_core::{
    causal::{self as core, Causal, GridVec, Resolution},
    time::{Duration, Timestamp},
    Error,
};

use crate::{
    convert::{err, invalid, type_name, IntoPy},
    time::{
        as_int, duration_from_py, period, uint_from_py, uint_to_py, Duration1024, Duration384,
        Duration512, Duration768, Given, Time1024, Time384, Time512, Time768, Width,
    },
};

/// A displacement argument: a sequence of exactly three `int`s, in grid steps.
pub(crate) struct Grid(GridVec);

impl<'a, 'py> FromPyObject<'a, 'py> for Grid {
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if obj.is_instance_of::<PyString>() || obj.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "expected a sequence of three ints (grid steps), not text",
            ));
        }
        let len = obj.len().map_err(|_| {
            PyTypeError::new_err(format!(
                "expected a sequence of three ints (grid steps), got {}",
                type_name(&obj)
            ))
        })?;
        if len != 3 {
            return Err(PyValueError::new_err(format!(
                "expected three components (x, y, z) in grid steps, got {len}"
            )));
        }
        // The extraction goes through `__index__`: an int or a NumPy integer passes, a float does not.
        let component = |i: usize| -> PyResult<i128> {
            let item = obj.get_item(i)?;
            item.extract::<i128>().map_err(|e| {
                if e.is_instance_of::<PyTypeError>(obj.py()) {
                    PyTypeError::new_err(format!(
                        "a displacement is whole grid steps as ints, never floats; got {}",
                        type_name(&item)
                    ))
                } else {
                    PyOverflowError::new_err(
                        "a displacement component must fit a signed 128-bit integer of grid steps",
                    )
                }
            })
        };
        Ok(Self(GridVec::new(
            component(0)?,
            component(1)?,
            component(2)?,
        )))
    }
}

/// The widths, in the order the error messages list them.
const WIDTHS: &str = "384, 512, 768 or 1024";

/// Call a function generic over the width with the limb count and the time class for `$bits`, which must be
/// one of the four widths.
macro_rules! at_width {
    ($bits:expr, $f:ident($($arg:expr),* $(,)?)) => {
        match $bits {
            384 => $f::<6, Time384>($($arg),*),
            512 => $f::<8, Time512>($($arg),*),
            768 => $f::<12, Time768>($($arg),*),
            1024 => $f::<16, Time1024>($($arg),*),
            other => unreachable!("no time class is {other} bits wide"),
        }
    };
}

/// The width in bits of a time object: the width a module-level function works at.
fn time_bits(obj: &Bound<'_, PyAny>) -> PyResult<u32> {
    if obj.is_instance_of::<Time384>() {
        Ok(384)
    } else if obj.is_instance_of::<Time512>() {
        Ok(512)
    } else if obj.is_instance_of::<Time768>() {
        Ok(768)
    } else if obj.is_instance_of::<Time1024>() {
        Ok(1024)
    } else {
        Err(PyTypeError::new_err(format!(
            "expected a Time384, Time512, Time768 or Time1024, got {}",
            type_name(obj)
        )))
    }
}

/// The width in bits of a duration object, or `None` for anything else.
fn duration_bits(obj: &Bound<'_, PyAny>) -> Option<u32> {
    if obj.is_instance_of::<Duration384>() {
        Some(384)
    } else if obj.is_instance_of::<Duration512>() {
        Some(512)
    } else if obj.is_instance_of::<Duration768>() {
        Some(768)
    } else if obj.is_instance_of::<Duration1024>() {
        Some(1024)
    } else {
        None
    }
}

/// A time argument of the width `W`. Times of different widths do not mix: any other object is a `TypeError`.
fn time_arg<const N: usize, W: Width<N>>(
    obj: &Bound<'_, PyAny>,
    what: &str,
) -> PyResult<Timestamp<N>> {
    W::time_from(obj).ok_or_else(|| {
        PyTypeError::new_err(format!(
            "expected a {} for {what} (times of different widths do not mix), got {}",
            W::TIME,
            type_name(obj)
        ))
    })
}

/// A duration argument of the width `W`: a duration of that width or an `int` of ticks.
fn duration_arg<const N: usize, W: Width<N>>(
    obj: &Bound<'_, PyAny>,
    what: &str,
) -> PyResult<Duration<N>> {
    if let Some(d) = W::duration_from(obj) {
        Ok(d)
    } else if duration_bits(obj).is_some() {
        Err(PyTypeError::new_err(format!(
            "expected a {} or an int of ticks for {what} (durations of different widths do not mix), got {}",
            W::DURATION,
            type_name(obj)
        )))
    } else {
        duration_from_py::<N>(obj, what)
    }
}

fn light_ticks_at<const N: usize, W: Width<N>>(
    py: Python<'_>,
    res: Resolution,
    dx: GridVec,
) -> PyResult<Bound<'_, PyAny>> {
    let ticks = core::Clock::<N>::new(res).light_ticks(dx).py()?;
    uint_to_py(py, ticks.magnitude_ticks())
}

fn first_arrival_at<'py, const N: usize, W: Width<N>>(
    res: Resolution,
    emission: &Bound<'py, PyAny>,
    dx: GridVec,
) -> PyResult<Bound<'py, PyAny>> {
    let e = time_arg::<N, W>(emission, "emission")?;
    let arrival = core::Clock::<N>::new(res).first_arrival(e, dx).py()?;
    W::time_to(emission.py(), arrival)
}

fn arrived_at<const N: usize, W: Width<N>>(
    res: Resolution,
    emission: &Bound<'_, PyAny>,
    now: &Bound<'_, PyAny>,
    dx: GridVec,
) -> PyResult<bool> {
    let e = time_arg::<N, W>(emission, "emission")?;
    let n = time_arg::<N, W>(now, "now")?;
    Ok(core::Clock::<N>::new(res).arrived(e, n, dx))
}

fn classify_at<const N: usize, W: Width<N>>(
    res: Resolution,
    dt: &Bound<'_, PyAny>,
    dx: GridVec,
) -> PyResult<&'static str> {
    let dt = duration_arg::<N, W>(dt, "dt")?;
    Ok(match core::Clock::<N>::new(res).classify(dt, dx) {
        Causal::Timelike => "timelike",
        Causal::Null => "null",
        Causal::Spacelike => "spacelike",
    })
}

/// Run the core's search with a Python callback, keeping the callback's own exception if it raises.
fn settle_arrival_at<'py, const N: usize, W: Width<N>>(
    res: Resolution,
    emission: &Bound<'py, PyAny>,
    estimate: &Bound<'py, PyAny>,
    displacement_at: &Bound<'py, PyAny>,
    guard: u32,
) -> PyResult<Bound<'py, PyAny>> {
    let py = emission.py();
    let e = time_arg::<N, W>(emission, "emission")?;
    let est = time_arg::<N, W>(estimate, "estimate")?;
    let displacement = |tick: Timestamp<N>| -> PyResult<GridVec> {
        let tick = W::time_to(py, tick)?;
        Ok(displacement_at.call1((tick,))?.extract::<Grid>()?.0)
    };
    let mut raised: Option<PyErr> = None;
    let result = core::Clock::<N>::new(res).settle_arrival(e, est, guard, |tick| {
        displacement(tick).map_err(|e| {
            raised = Some(e);
            Error::InvalidInput
        })
    });
    match (result, raised) {
        (_, Some(e)) => Err(e),
        (Ok(tick), None) => W::time_to(py, tick),
        (Err(e), None) => Err(err(e)),
    }
}

fn time_of<'py, const N: usize, W: Width<N>>(
    ticks: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let value = uint_from_py::<N>(ticks, "ticks")?;
    W::time_to(ticks.py(), Timestamp::from_raw_ticks(value))
}

fn duration_of<'py, const N: usize, W: Width<N>>(
    ticks: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let value = duration_from_py::<N>(ticks, "ticks")?;
    W::duration_to(ticks.py(), value)
}

/// A whole number of seconds as a duration. The core's `Clock::seconds` takes a `u64`, which stops at 584
/// billion years; this multiplies at the full width and keeps the sign.
fn seconds_at<'py, const N: usize, W: Width<N>>(
    res: Resolution,
    seconds: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let seconds_given = duration_from_py::<N>(seconds, "seconds")?;
    let ticks = seconds_given
        .magnitude_ticks()
        .checked_mul(BUint::from_digit(res.ticks_per_second()))
        .ok_or(Error::Overflow)
        .py()?;
    let duration = Duration::from_parts(seconds_given.is_negative(), ticks);
    W::duration_to(seconds.py(), duration)
}

fn local_seconds_at<const N: usize, W: Width<N>>(
    res: Resolution,
    duration: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    let duration = duration_arg::<N, W>(duration, "duration")?;
    core::Clock::<N>::new(res).local_seconds(duration).py()
}

fn local_seconds_between_at<const N: usize, W: Width<N>>(
    res: Resolution,
    later: &Bound<'_, PyAny>,
    earlier: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    let later = time_arg::<N, W>(later, "later")?;
    let earlier = time_arg::<N, W>(earlier, "earlier")?;
    core::Clock::<N>::new(res)
        .local_seconds_between(later, earlier)
        .py()
}

fn phase_seconds_at<const N: usize, W: Width<N>>(
    res: Resolution,
    t: &Bound<'_, PyAny>,
    period_ticks: u64,
) -> PyResult<f64> {
    let t = time_arg::<N, W>(t, "t")?;
    core::Clock::<N>::new(res)
        .phase_seconds(t, period_ticks)
        .py()
}

fn reach_at<const N: usize, W: Width<N>>(res: Resolution) -> f64 {
    core::Clock::<N>::new(res).reach_log10_years()
}

fn time_class<const N: usize, W: Width<N>>(py: Python<'_>) -> Bound<'_, PyType> {
    W::time_class(py)
}

fn duration_class<const N: usize, W: Width<N>>(py: Python<'_>) -> Bound<'_, PyType> {
    W::duration_class(py)
}

/// A whole-number argument of the `Clock` constructor: `default` if it was left out, `None` if it is negative
/// or beyond 64 bits. Anything but an `int` is a `TypeError`.
fn count(given: &Given<'_>, name: &str, default: u64) -> PyResult<Option<u64>> {
    let Some(obj) = &given.0 else {
        return Ok(Some(default));
    };
    let int = as_int(obj).ok_or_else(|| {
        PyTypeError::new_err(format!("{name} must be an int, got {}", type_name(obj)))
    })?;
    Ok(int.extract::<u64>().ok())
}

/// A constructor argument as the error messages show it.
fn shown(given: &Given<'_>) -> String {
    given
        .0
        .as_ref()
        .and_then(|obj| obj.repr().ok())
        .map_or_else(|| "?".to_owned(), |repr| repr.to_string())
}

/// A stored clock: a width, a tick and a grid spacing, chosen together.
///
/// `Clock(bits=384, ticks_per_second=1000, grid_per_metre=1)`. `bits` is the width of the stored integer: 384,
/// 512, 768 or 1024. `ticks_per_second` is how many ticks make a second, from 1 to 10**12 (a picosecond).
/// `grid_per_metre` is how many grid steps make a metre, from 1 to 10**6 (a micrometre). The default is a
/// millisecond tick on a metre grid, which is what the module-level functions of `srkit.causal` use.
///
/// The times of a clock are the time class of its width (`clock.Time`, built with `clock.time(ticks)`), and its
/// durations are the matching duration class. A time counts ticks and does not carry the clock, so keep one
/// clock per application and use it for every call. Times of another width are refused with `TypeError`.
///
/// Every decision is exact. Displacements are three ints of grid steps (`clock.grid` rounds metres to them),
/// and with `f` ticks per second and `g` grid steps per metre
///
///     arrived  <=>  c * g * dt_ticks >= f * |dx_steps|
///
/// holds or fails with no rounding. The grid should be much finer than the distance light covers in a tick
/// (`light_steps_per_tick` in the hundreds or more), or rounding a position to the grid moves an arrival by
/// whole ticks: a nanosecond tick wants a millimetre grid, not a metre one.
///
/// Clocks compare equal, and hash alike, when all three numbers match.
#[pyclass(frozen, module = "srkit.causal", name = "Clock")]
#[derive(Clone, Copy)]
pub struct Clock {
    bits: u32,
    res: Resolution,
}

impl Clock {
    fn key(&self) -> (u32, u64, u64) {
        (
            self.bits,
            self.res.ticks_per_second(),
            self.res.grid_per_metre(),
        )
    }
}

#[pymethods]
impl Clock {
    #[new]
    #[pyo3(
        signature = (bits=Given(None), ticks_per_second=Given(None), grid_per_metre=Given(None)),
        text_signature = "(bits=384, ticks_per_second=1000, grid_per_metre=1)"
    )]
    fn new(
        bits: Given<'_>,
        ticks_per_second: Given<'_>,
        grid_per_metre: Given<'_>,
    ) -> PyResult<Self> {
        let width = match count(&bits, "bits", 384)? {
            Some(width @ (384 | 512 | 768 | 1024)) => width as u32,
            _ => {
                return Err(invalid(format!(
                    "bits must be {WIDTHS}, got {}",
                    shown(&bits)
                )))
            }
        };
        let ticks = count(&ticks_per_second, "ticks_per_second", 1000)?;
        let grid = count(&grid_per_metre, "grid_per_metre", 1)?;
        let res = match (ticks, grid) {
            (Some(ticks), Some(grid)) => Resolution::new(ticks, grid).ok(),
            _ => None,
        };
        let res = res.ok_or_else(|| {
            let or_default = |given: &Given<'_>, default: &str| {
                if given.0.is_some() {
                    shown(given)
                } else {
                    default.to_owned()
                }
            };
            invalid(format!(
                "ticks_per_second must be from 1 to 10**12 and grid_per_metre from 1 to 10**6, got {} and {}",
                or_default(&ticks_per_second, "1000"),
                or_default(&grid_per_metre, "1"),
            ))
        })?;
        Ok(Self { bits: width, res })
    }

    /// Width of the stored integer in bits: 384, 512, 768 or 1024.
    #[getter]
    fn bits(&self) -> u32 {
        self.bits
    }
    /// How many ticks make one second.
    #[getter]
    fn ticks_per_second(&self) -> u64 {
        self.res.ticks_per_second()
    }
    /// How many grid steps make one metre.
    #[getter]
    fn grid_per_metre(&self) -> u64 {
        self.res.grid_per_metre()
    }
    /// How many grid steps light covers in one tick. Rounding a position to the grid moves an arrival by up to
    /// one part in this many of a tick, so it should be large (hundreds or more).
    #[getter]
    fn light_steps_per_tick(&self) -> f64 {
        self.res.light_steps_per_tick()
    }
    /// The length of one tick in seconds.
    #[getter]
    fn tick_seconds(&self) -> f64 {
        self.res.tick_seconds()
    }
    /// The length of one grid step in metres.
    #[getter]
    fn grid_metres(&self) -> f64 {
        self.res.grid_metres()
    }
    /// The base-10 logarithm of the longest time the clock can hold, in Julian years: 100 means 10**100 years.
    #[getter]
    fn reach_log10_years(&self) -> f64 {
        at_width!(self.bits, reach_at(self.res))
    }
    /// The time class of this width: `srkit.Time384`, `Time512`, `Time768` or `Time1024`.
    #[getter]
    #[pyo3(name = "Time")]
    fn time_type<'py>(&self, py: Python<'py>) -> Bound<'py, PyType> {
        at_width!(self.bits, time_class(py))
    }
    /// The duration class of this width: `srkit.Duration384`, `Duration512`, `Duration768` or `Duration1024`.
    #[getter]
    #[pyo3(name = "Duration")]
    fn duration_type<'py>(&self, py: Python<'py>) -> Bound<'py, PyType> {
        at_width!(self.bits, duration_class(py))
    }

    /// The time `ticks` ticks after the epoch, as the time class of this width. `ticks` is an `int` that is
    /// not negative. Raises `OverflowError` if it does not fit the width.
    fn time<'py>(&self, ticks: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        at_width!(self.bits, time_of(ticks))
    }
    /// A duration of `ticks` ticks, as the duration class of this width. `ticks` is an `int` of either sign.
    fn duration<'py>(&self, ticks: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        at_width!(self.bits, duration_of(ticks))
    }
    /// A duration of `n` whole seconds: `n * ticks_per_second` ticks. `n` is an `int` of either sign and any
    /// size the width can hold.
    fn seconds<'py>(&self, n: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        at_width!(self.bits, seconds_at(self.res, n))
    }
    /// Round a displacement in metres, given as floats, to this clock's grid (to nearest, ties away from
    /// zero). Returns three ints of grid steps.
    ///
    /// A float holds 53 bits, so beyond 9e15 grid steps it no longer names a single step: keep far positions
    /// as ints. Raises `OverflowError` for a component that is not finite or is beyond the grid.
    fn grid(&self, x_m: f64, y_m: f64, z_m: f64) -> PyResult<(i128, i128, i128)> {
        let g = GridVec::from_metres_at(self.res, x_m, y_m, z_m).py()?;
        Ok((g.x, g.y, g.z))
    }

    /// The light-travel time over the displacement `dx` (three ints, grid steps), rounded up to whole ticks.
    /// Exact, by an integer square root. Returns an `int`.
    fn light_ticks<'py>(&self, py: Python<'py>, dx: Grid) -> PyResult<Bound<'py, PyAny>> {
        at_width!(self.bits, light_ticks_at(py, self.res, dx.0))
    }
    /// The first tick at which a signal emitted at `emission` has reached a fixed point `dx` away (three ints,
    /// grid steps). Returns a time of this clock's width.
    fn first_arrival<'py>(
        &self,
        emission: &Bound<'py, PyAny>,
        dx: Grid,
    ) -> PyResult<Bound<'py, PyAny>> {
        at_width!(self.bits, first_arrival_at(self.res, emission, dx.0))
    }
    /// Whether a signal emitted at `emission` has reached a point `dx` away (three ints, grid steps) by the
    /// tick `now`. Exact: true exactly from `first_arrival(emission, dx)` onwards, and false before the
    /// emission.
    fn arrived(
        &self,
        emission: &Bound<'_, PyAny>,
        now: &Bound<'_, PyAny>,
        dx: Grid,
    ) -> PyResult<bool> {
        at_width!(self.bits, arrived_at(self.res, emission, now, dx.0))
    }
    /// Classify the separation between two events `dt` apart in time (an `int` of ticks or a duration of this
    /// clock's width; the sign does not matter) and `dx` apart in space (three ints, grid steps).
    ///
    /// Returns 'timelike', 'null' or 'spacelike'. Exact.
    fn classify(&self, dt: &Bound<'_, PyAny>, dx: Grid) -> PyResult<&'static str> {
        at_width!(self.bits, classify_at(self.res, dt, dx.0))
    }
    /// The first tick at which a signal reaches a moving receiver, settled exactly from a floating-point
    /// estimate.
    ///
    /// `displacement_at(tick)` is called with a time of this clock's width and must return the receiver's
    /// position relative to the emission point at that tick, as three ints of grid steps (evaluate the
    /// receiver's worldline in floats and round with `grid`). Starting from `estimate`, the search steps back
    /// while the signal has already arrived and forward until it has, so the result is the first tick that
    /// passes the exact check and can never be a tick early. `guard` bounds the number of steps each way; an
    /// estimate further off than that raises `srkit.PrecisionLoss`.
    #[pyo3(signature = (emission, estimate, displacement_at, guard=8))]
    fn settle_arrival<'py>(
        &self,
        emission: &Bound<'py, PyAny>,
        estimate: &Bound<'py, PyAny>,
        displacement_at: &Bound<'py, PyAny>,
        guard: u32,
    ) -> PyResult<Bound<'py, PyAny>> {
        at_width!(
            self.bits,
            settle_arrival_at(self.res, emission, estimate, displacement_at, guard)
        )
    }

    /// A duration (of this clock's width, or an `int` of ticks) in seconds, as a float. Raises
    /// `srkit.PrecisionLoss` beyond 2**53 ticks, where a float can no longer hold every tick.
    fn local_seconds(&self, duration: &Bound<'_, PyAny>) -> PyResult<f64> {
        at_width!(self.bits, local_seconds_at(self.res, duration))
    }
    /// `later - earlier` in seconds. The subtraction is exact; the result is a float. Raises
    /// `srkit.PrecisionLoss` if the difference is beyond 2**53 ticks.
    fn local_seconds_between(
        &self,
        later: &Bound<'_, PyAny>,
        earlier: &Bound<'_, PyAny>,
    ) -> PyResult<f64> {
        at_width!(
            self.bits,
            local_seconds_between_at(self.res, later, earlier)
        )
    }
    /// The exact phase of the time `t` in a period of `period_ticks` whole ticks, in seconds. The reduction is
    /// in integers, so it is right even at huge epochs.
    fn phase_seconds(&self, t: &Bound<'_, PyAny>, period_ticks: Bound<'_, PyInt>) -> PyResult<f64> {
        let period_ticks = period(&period_ticks)?;
        at_width!(self.bits, phase_seconds_at(self.res, t, period_ticks))
    }

    fn __repr__(&self) -> String {
        let (bits, ticks, grid) = self.key();
        format!("Clock(bits={bits}, ticks_per_second={ticks}, grid_per_metre={grid})")
    }
    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        self.key().into_pyobject(py)?.hash()
    }
    fn __richcmp__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
        op: CompareOp,
    ) -> PyResult<Bound<'py, PyAny>> {
        let same = match other.cast::<Clock>() {
            Ok(other) => self.key() == other.get().key(),
            Err(_) => return Ok(py.NotImplemented().into_bound(py)),
        };
        let result = match op {
            CompareOp::Eq => same,
            CompareOp::Ne => !same,
            _ => return Ok(py.NotImplemented().into_bound(py)),
        };
        Ok(result.into_pyobject(py)?.to_owned().into_any())
    }
}

/// The resolution of the module-level functions: a millisecond tick on a metre grid.
const DEFAULT: Resolution = Resolution::MILLISECOND_METRE;

/// The light-travel time over the displacement `dx` (three ints, metres), rounded up to whole milliseconds:
/// the smallest `n` with `c n >= 1000 |dx|`. Exact, by an integer square root. Returns an `int`.
///
/// Exactly one light-second (299,792,458 m) is 1000 ticks; one metre more is 1001. For another tick or grid
/// use `Clock.light_ticks`.
#[pyfunction]
fn light_ticks(py: Python<'_>, dx: Grid) -> PyResult<Bound<'_, PyAny>> {
    light_ticks_at::<8, Time512>(py, DEFAULT, dx.0)
}

/// The first tick at which a signal emitted at `emission` has reached a fixed point `dx` away (three ints,
/// metres), with millisecond ticks. Returns a time of the same width as `emission`.
///
/// The same as `Clock(bits).first_arrival`, where `bits` is the width of `emission`.
#[pyfunction]
fn first_arrival<'py>(emission: &Bound<'py, PyAny>, dx: Grid) -> PyResult<Bound<'py, PyAny>> {
    at_width!(
        time_bits(emission)?,
        first_arrival_at(DEFAULT, emission, dx.0)
    )
}

/// Whether a signal emitted at `emission` has reached a point `dx` away (three ints, metres) by the tick
/// `now`, with millisecond ticks. Exact: true exactly from `first_arrival(emission, dx)` onwards, and false
/// before the emission.
///
/// The same as `Clock(bits).arrived`, where `bits` is the width of the two times.
#[pyfunction]
fn arrived(emission: &Bound<'_, PyAny>, now: &Bound<'_, PyAny>, dx: Grid) -> PyResult<bool> {
    at_width!(
        time_bits(emission)?,
        arrived_at(DEFAULT, emission, now, dx.0)
    )
}

/// Classify the separation between two events `dt` apart in time (an `int` of milliseconds or a duration of
/// any width; the sign does not matter) and `dx` apart in space (three ints, metres).
///
/// Returns 'timelike' (the later event is strictly inside the earlier one's light cone), 'null' (light from
/// the earlier reaches the later exactly) or 'spacelike' (neither can influence the other). Exact.
///
/// The same as `Clock(bits).classify`, where `bits` is the width of the duration; an `int` is read at 512
/// bits.
#[pyfunction]
fn classify(dt: &Bound<'_, PyAny>, dx: Grid) -> PyResult<&'static str> {
    at_width!(
        duration_bits(dt).unwrap_or(512),
        classify_at(DEFAULT, dt, dx.0)
    )
}

/// Round a displacement in metres, given as floats, to the one-metre grid (to nearest, ties away from zero).
/// Returns three ints. Raises `OverflowError` for a component that is not finite or is beyond the grid.
///
/// For a grid of another spacing use `Clock.grid`.
#[pyfunction]
fn grid_from_metres(x: f64, y: f64, z: f64) -> PyResult<(i128, i128, i128)> {
    let g = GridVec::from_metres(x, y, z).py()?;
    Ok((g.x, g.y, g.z))
}

/// The first tick at which a signal reaches a moving receiver, settled exactly from a floating-point estimate,
/// with millisecond ticks on a metre grid.
///
/// `displacement_at(tick)` is called with a time and must return the receiver's position relative to the
/// emission point at that tick, as three ints in metres (evaluate the receiver's worldline in floats and round
/// with `grid_from_metres`). Starting from `estimate`, the search steps back while the signal has already
/// arrived and forward until it has, so the result is the first tick that passes the exact check and can never
/// be a tick early. `guard` bounds the number of steps each way; an estimate further off than that raises
/// `srkit.PrecisionLoss`. Returns a time of the same width as `emission`.
///
/// The same as `Clock(bits).settle_arrival`, where `bits` is the width of the two times.
#[pyfunction]
#[pyo3(signature = (emission, estimate, displacement_at, guard=8))]
fn settle_arrival<'py>(
    emission: &Bound<'py, PyAny>,
    estimate: &Bound<'py, PyAny>,
    displacement_at: &Bound<'py, PyAny>,
    guard: u32,
) -> PyResult<Bound<'py, PyAny>> {
    at_width!(
        time_bits(emission)?,
        settle_arrival_at(DEFAULT, emission, estimate, displacement_at, guard)
    )
}

const DOC: &str = "\
Exact causal decisions on stored integer times.

Geometry in srkit is floating point, which is right for drawing and wrong for deciding. \"Has this signal
arrived by tick n?\" and \"is event B inside A's light cone?\" are answered here in integer arithmetic: times
are whole ticks (`srkit.Time384`, `Time512`, `Time768`, `Time1024`) and displacements are whole grid steps,
given as three Python ints (each within the range of a signed 128-bit integer). The speed of light is exactly
299,792,458 m/s, so with millisecond ticks on a metre grid

    arrived  <=>  c * dt_ms >= 1000 * |dx_m|

holds or fails with no rounding at all. Compute where things are in floating point, round the endpoints to the
grid, then decide exactly. Floats are refused as displacements.

The tick and the grid are the application's choice. `Clock(bits, ticks_per_second, grid_per_metre)` holds a
width (384, 512, 768 or 1024 bits), a tick and a grid spacing, and its methods are the checks at that
resolution. The functions of this module are the default, a millisecond tick on a metre grid: each is the
method of the same name on `Clock(bits, 1000, 1)`, where `bits` is the width of the times it is given, and
`grid_from_metres` is `Clock().grid`. On the metre grid a displacement spans 1.8e22 light-years each way.

Three resolutions are named, as `(ticks_per_second, grid_per_metre)` pairs to unpack into `Clock`:
`MILLISECOND_METRE` (1000, 1), `MICROSECOND_METRE` (10**6, 1) and `NANOSECOND_MILLIMETRE` (10**9, 1000), as
in `Clock(768, *NANOSECOND_MILLIMETRE)`. `C_METRES_PER_SECOND` is the speed of light, an int.
";

/// A named resolution as the pair `Clock` takes.
fn pair(res: Resolution) -> (u64, u64) {
    (res.ticks_per_second(), res.grid_per_metre())
}

pub(crate) fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.causal")?;
    m.setattr("__doc__", DOC)?;
    m.add("C_METRES_PER_SECOND", core::C_METRES_PER_SECOND)?;
    m.add("MILLISECOND_METRE", pair(Resolution::MILLISECOND_METRE))?;
    m.add("MICROSECOND_METRE", pair(Resolution::MICROSECOND_METRE))?;
    m.add(
        "NANOSECOND_MILLIMETRE",
        pair(Resolution::NANOSECOND_MILLIMETRE),
    )?;
    m.add_class::<Clock>()?;
    m.add_function(wrap_pyfunction!(light_ticks, &m)?)?;
    m.add_function(wrap_pyfunction!(first_arrival, &m)?)?;
    m.add_function(wrap_pyfunction!(arrived, &m)?)?;
    m.add_function(wrap_pyfunction!(classify, &m)?)?;
    m.add_function(wrap_pyfunction!(grid_from_metres, &m)?)?;
    m.add_function(wrap_pyfunction!(settle_arrival, &m)?)?;
    m.setattr(
        "__all__",
        [
            "C_METRES_PER_SECOND",
            "Clock",
            "MICROSECOND_METRE",
            "MILLISECOND_METRE",
            "NANOSECOND_MILLIMETRE",
            "arrived",
            "classify",
            "first_arrival",
            "grid_from_metres",
            "light_ticks",
            "settle_arrival",
        ],
    )?;
    Ok(m)
}
