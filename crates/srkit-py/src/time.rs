//! `srkit.Time384`, `Time512`, `Time768`, `Time1024` and their durations: exact integer ticks, to and from
//! Python `int`.
//!
//! A time counts whole ticks and does not know how long a tick is. A tick is a millisecond by default, and the
//! methods on these classes that give seconds read it so; `srkit.causal.Clock` (in `causal.rs`) holds a width
//! together with a tick and a grid spacing, and reaches the classes through [`Width`].
//!
//! A Python `int` has no fixed width, so the conversion is lossless both ways: through `u64` when the value is
//! small, through little-endian bytes otherwise. Nothing here goes through a float.

use bnum::BUint;
use pyo3::{
    basic::CompareOp,
    exceptions::{PyOverflowError, PyTypeError},
    prelude::*,
    types::{PyBytes, PyInt, PyType},
};
use srkit_core::time::{Duration, Timestamp};

use crate::convert::{invalid, type_name, IntoPy};

/// A non-negative Python `int` as a fixed-width unsigned integer.
///
/// `TypeError` for anything but an `int`, `ValueError` for a negative one, `OverflowError` if it needs more
/// than `64 N` bits.
pub(crate) fn uint_from_py<const N: usize>(
    obj: &Bound<'_, PyAny>,
    what: &str,
) -> PyResult<BUint<N>> {
    let int = int_from_py(obj, what)?;
    if int.lt(0)? {
        return Err(invalid(format!("{what} must not be negative")));
    }
    magnitude_from_py(int.as_any(), what)
}

/// An exact integer argument as a Python `int`: an `int`, or an integer type that says so through `__index__`
/// (a NumPy integer, say). `None` for anything else, floats included.
pub(crate) fn as_int<'py>(obj: &Bound<'py, PyAny>) -> Option<Bound<'py, PyInt>> {
    if let Ok(int) = obj.cast::<PyInt>() {
        return Some(int.clone());
    }
    obj.call_method0("__index__")
        .ok()?
        .cast_into::<PyInt>()
        .ok()
}

fn int_from_py<'py>(obj: &Bound<'py, PyAny>, what: &str) -> PyResult<Bound<'py, PyInt>> {
    as_int(obj).ok_or_else(|| {
        PyTypeError::new_err(format!(
            "{what} must be an int (an exact count, never a float), got {}",
            type_name(obj)
        ))
    })
}

/// A non-negative Python `int`, already checked, as `N` 64-bit digits.
fn magnitude_from_py<const N: usize>(int: &Bound<'_, PyAny>, what: &str) -> PyResult<BUint<N>> {
    if let Ok(small) = int.extract::<u64>() {
        return Ok(BUint::from_digit(small));
    }
    let bits: usize = int.call_method0("bit_length")?.extract()?;
    if bits > N * 64 {
        return Err(PyOverflowError::new_err(format!(
            "{what} needs {bits} bits and does not fit in {}",
            N * 64
        )));
    }
    let bytes = int.call_method1("to_bytes", (N * 8, "little"))?;
    let bytes = bytes.cast::<PyBytes>()?.as_bytes();
    let mut digits = [0u64; N];
    for (digit, chunk) in digits.iter_mut().zip(bytes.chunks_exact(8)) {
        *digit = u64::from_le_bytes(chunk.try_into().expect("chunks of eight bytes"));
    }
    Ok(BUint::from_digits(digits))
}

/// A fixed-width unsigned integer as a Python `int`.
pub(crate) fn uint_to_py<const N: usize>(
    py: Python<'_>,
    value: BUint<N>,
) -> PyResult<Bound<'_, PyAny>> {
    let digits = value.digits();
    if digits[1..].iter().all(|d| *d == 0) {
        return Ok(digits[0].into_pyobject(py)?.into_any());
    }
    let mut bytes = Vec::with_capacity(N * 8);
    for digit in digits {
        bytes.extend_from_slice(&digit.to_le_bytes());
    }
    py.get_type::<PyInt>()
        .call_method1("from_bytes", (PyBytes::new(py, &bytes), "little"))
}

/// A Python `int` of either sign as a signed duration.
pub(crate) fn duration_from_py<const N: usize>(
    obj: &Bound<'_, PyAny>,
    what: &str,
) -> PyResult<Duration<N>> {
    let int = int_from_py(obj, what)?;
    let negative = int.lt(0)?;
    let magnitude = if negative {
        magnitude_from_py::<N>(&int.neg()?, what)?
    } else {
        magnitude_from_py::<N>(int.as_any(), what)?
    };
    Ok(Duration::from_parts(negative, magnitude))
}

/// A signed duration as a Python `int`.
pub(crate) fn duration_to_py<const N: usize>(
    py: Python<'_>,
    value: Duration<N>,
) -> PyResult<Bound<'_, PyAny>> {
    let magnitude = uint_to_py(py, value.magnitude_millis())?;
    if value.is_negative() {
        magnitude.neg()
    } else {
        Ok(magnitude)
    }
}

/// Order two signed durations.
fn compare<const N: usize>(a: Duration<N>, b: Duration<N>) -> std::cmp::Ordering {
    match (a.is_negative(), b.is_negative()) {
        (false, false) => a.magnitude_millis().cmp(&b.magnitude_millis()),
        (true, true) => b.magnitude_millis().cmp(&a.magnitude_millis()),
        (false, true) => std::cmp::Ordering::Greater,
        (true, false) => std::cmp::Ordering::Less,
    }
}

/// A constructor argument that may be left out. Unlike `Option`, an explicit `None` is not "left out": it is an
/// object like any other, and is refused as not being an `int`.
pub(crate) struct Given<'py>(pub(crate) Option<Bound<'py, PyAny>>);

impl<'a, 'py> FromPyObject<'a, 'py> for Given<'py> {
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        Ok(Self(Some(obj.to_owned())))
    }
}

/// A period in whole ticks, for modular phase.
pub(crate) fn period(ticks: &Bound<'_, PyInt>) -> PyResult<u64> {
    if ticks.le(0)? {
        return Err(invalid("the period must be a positive number of ticks"));
    }
    ticks
        .extract::<u64>()
        .map_err(|_| PyOverflowError::new_err("the period must be below 2**64 ticks"))
}

/// One width of exact time: the Python classes that hold `Timestamp<N>` and `Duration<N>`.
///
/// Implemented by each time class, so that code generic over the width (`srkit.causal`) can read and build
/// times and durations without naming the classes.
pub(crate) trait Width<const N: usize> {
    /// The Python name of the time class.
    const TIME: &'static str;
    /// The Python name of the duration class.
    const DURATION: &'static str;
    /// The time an object holds, if it is a time of this width.
    fn time_from(obj: &Bound<'_, PyAny>) -> Option<Timestamp<N>>;
    /// The duration an object holds, if it is a duration of this width.
    fn duration_from(obj: &Bound<'_, PyAny>) -> Option<Duration<N>>;
    /// A new time object of this width.
    fn time_to(py: Python<'_>, value: Timestamp<N>) -> PyResult<Bound<'_, PyAny>>;
    /// A new duration object of this width.
    fn duration_to(py: Python<'_>, value: Duration<N>) -> PyResult<Bound<'_, PyAny>>;
    /// The time class.
    fn time_class(py: Python<'_>) -> Bound<'_, PyType>;
    /// The duration class.
    fn duration_class(py: Python<'_>) -> Bound<'_, PyType>;
}

macro_rules! exact_time {
    ($time:ident, $dur:ident, $n:literal, $tname:literal, $dname:literal) => {
        /// An exact time: whole ticks since an epoch of the application's choosing, held as an unsigned
        /// integer of fixed width (`BITS`: 384, 512, 768 or 1024). It converts to and from a Python `int`
        /// without loss.
        ///
        /// A time counts ticks and does not know how long a tick is. A tick is a millisecond by default: the
        /// methods here that give seconds or take metres (`local_seconds_since`, `phase_seconds`,
        /// `stationary_arrival`, `stationary_signal_arrived`) read it so. A `srkit.causal.Clock` gives ticks
        /// another length and does the same work at that length.
        ///
        /// Epochs are never turned into floats. Subtract two times first (the difference is an exact
        /// duration), then ask for `local_seconds_since` when the difference is small enough to be a float.
        ///
        /// Arithmetic: `time + d` and `time - d` take a duration or an `int` of ticks and give a time;
        /// `time - time` gives a duration; `time % period` is the exact phase in a period of whole ticks.
        /// Times of one width compare and hash by value, and times of different widths do not mix. A result
        /// before the epoch or past `MAX` raises `OverflowError`.
        #[pyclass(frozen, module = "srkit", name = $tname)]
        #[derive(Clone, Copy)]
        pub struct $time {
            pub(crate) inner: Timestamp<$n>,
        }

        /// An exact signed duration: whole ticks, with a magnitude as wide as the matching time type (`BITS`:
        /// 384, 512, 768 or 1024). It converts to and from a Python `int` without loss, and compares with
        /// durations of the same width and with `int`.
        ///
        /// A tick is a millisecond by default, which is how `local_seconds` and `light_travel_time` read it.
        /// A `srkit.causal.Clock` gives ticks another length.
        #[pyclass(frozen, module = "srkit", name = $dname)]
        #[derive(Clone, Copy)]
        pub struct $dur {
            pub(crate) inner: Duration<$n>,
        }

        impl $dur {
            /// A duration argument: a duration of this width or an `int` of ticks.
            fn arg(obj: &Bound<'_, PyAny>) -> PyResult<Option<Duration<$n>>> {
                if let Ok(d) = obj.cast::<$dur>() {
                    Ok(Some(d.get().inner))
                } else if let Some(int) = as_int(obj) {
                    duration_from_py::<$n>(int.as_any(), "the duration").map(Some)
                } else {
                    Ok(None)
                }
            }
        }

        impl Width<$n> for $time {
            const TIME: &'static str = $tname;
            const DURATION: &'static str = $dname;
            fn time_from(obj: &Bound<'_, PyAny>) -> Option<Timestamp<$n>> {
                obj.cast::<$time>().ok().map(|t| t.get().inner)
            }
            fn duration_from(obj: &Bound<'_, PyAny>) -> Option<Duration<$n>> {
                obj.cast::<$dur>().ok().map(|d| d.get().inner)
            }
            fn time_to(py: Python<'_>, value: Timestamp<$n>) -> PyResult<Bound<'_, PyAny>> {
                Ok(Bound::new(py, $time { inner: value })?.into_any())
            }
            fn duration_to(py: Python<'_>, value: Duration<$n>) -> PyResult<Bound<'_, PyAny>> {
                Ok(Bound::new(py, $dur { inner: value })?.into_any())
            }
            fn time_class(py: Python<'_>) -> Bound<'_, PyType> {
                py.get_type::<$time>()
            }
            fn duration_class(py: Python<'_>) -> Bound<'_, PyType> {
                py.get_type::<$dur>()
            }
        }

        #[pymethods]
        impl $time {
            #[new]
            #[pyo3(signature = (millis=Given(None)), text_signature = "(millis=0)")]
            fn new(millis: Given<'_>) -> PyResult<Self> {
                let inner = match &millis.0 {
                    None => Timestamp::ZERO,
                    Some(obj) => match obj.cast::<$time>() {
                        Ok(same) => same.get().inner,
                        Err(_) => Timestamp::from_raw_millis(uint_from_py::<$n>(obj, "millis")?),
                    },
                };
                Ok(Self { inner })
            }

            /// Width of the stored integer in bits.
            #[classattr]
            const BITS: u32 = $n * 64;

            /// The epoch: zero ticks.
            #[classattr]
            #[pyo3(name = "ZERO")]
            fn zero() -> Self {
                Self {
                    inner: Timestamp::ZERO,
                }
            }
            /// The largest time this width can hold.
            #[classattr]
            #[pyo3(name = "MAX")]
            fn max() -> Self {
                Self {
                    inner: Timestamp::MAX,
                }
            }

            /// Parse unsigned decimal digits. Signs, spaces and separators are rejected.
            #[staticmethod]
            fn from_decimal(text: &str) -> PyResult<Self> {
                Ok(Self {
                    inner: Timestamp::from_decimal(text).py()?,
                })
            }

            fn __int__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
                uint_to_py(py, self.inner.raw_millis())
            }
            fn __str__(&self) -> String {
                self.inner.to_string()
            }
            fn __repr__(&self) -> String {
                format!("{}({})", $tname, self.inner)
            }
            fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
                uint_to_py(py, self.inner.raw_millis())?.hash()
            }
            fn __richcmp__(&self, other: PyRef<'_, Self>, op: CompareOp) -> bool {
                op.matches(self.inner.cmp(&other.inner))
            }

            fn __add__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                match $dur::arg(other)? {
                    Some(d) => {
                        let inner = self.inner.checked_add(d).py()?;
                        Ok(Bound::new(py, Self { inner })?.into_any())
                    }
                    None => Ok(py.NotImplemented().into_bound(py)),
                }
            }
            fn __radd__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                self.__add__(py, other)
            }
            fn __sub__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                if let Ok(origin) = other.cast::<$time>() {
                    let d = self.inner.signed_duration_since(origin.get().inner);
                    return Ok(Bound::new(py, $dur { inner: d })?.into_any());
                }
                match $dur::arg(other)? {
                    Some(d) => {
                        let inner = self.inner.checked_sub(d).py()?;
                        Ok(Bound::new(py, Self { inner })?.into_any())
                    }
                    None => Ok(py.NotImplemented().into_bound(py)),
                }
            }
            fn __mod__(&self, period_ms: Bound<'_, PyInt>) -> PyResult<u64> {
                self.inner.modulo_millis(period(&period_ms)?).py()
            }

            /// `self - origin` in seconds, reading a tick as a millisecond. The subtraction is exact; the
            /// result is a float. For ticks of another length use `Clock.local_seconds_between`.
            ///
            /// Raises `srkit.PrecisionLoss` if the difference is beyond 2**53 ticks (285,000 years of
            /// milliseconds), where a float can no longer hold every tick.
            fn local_seconds_since(&self, origin: PyRef<'_, Self>) -> PyResult<f64> {
                self.inner.local_seconds_since(origin.inner).py()
            }
            /// The exact phase in a period of `period_ms` whole ticks, in seconds, reading a tick as a
            /// millisecond. The reduction is in integers, so it is right even at huge epochs. For ticks of
            /// another length use `Clock.phase_seconds`.
            fn phase_seconds(&self, period_ms: Bound<'_, PyInt>) -> PyResult<f64> {
                self.inner.phase_seconds(period(&period_ms)?).py()
            }
            /// `self - origin` as a duration that is not negative. Raises `ValueError` if `origin` is later.
            fn duration_since(&self, origin: PyRef<'_, Self>) -> PyResult<$dur> {
                Ok($dur {
                    inner: self.inner.duration_since(origin.inner).py()?,
                })
            }
            /// The first tick at which a signal emitted now has crossed `distance_m` whole metres to a fixed
            /// point: the light-travel time rounded up to a whole tick, reading a tick as a millisecond.
            /// Exact. For another tick or grid use `Clock.first_arrival`.
            fn stationary_arrival(&self, distance_m: &Bound<'_, PyAny>) -> PyResult<Self> {
                let distance = uint_from_py::<$n>(distance_m, "distance_m")?;
                Ok(Self {
                    inner: self.inner.stationary_arrival_metres(distance).py()?,
                })
            }
            /// Whether this tick is late enough for a signal emitted at `emission` to have crossed
            /// `distance_m` whole metres to a fixed point, reading a tick as a millisecond. Exact; false
            /// before the emission. For another tick or grid use `Clock.arrived`.
            fn stationary_signal_arrived(
                &self,
                emission: PyRef<'_, Self>,
                distance_m: &Bound<'_, PyAny>,
            ) -> PyResult<bool> {
                let distance = uint_from_py::<$n>(distance_m, "distance_m")?;
                self.inner
                    .stationary_signal_arrived(emission.inner, distance)
                    .py()
            }
        }

        #[pymethods]
        impl $dur {
            #[new]
            #[pyo3(signature = (millis=Given(None)), text_signature = "(millis=0)")]
            fn new(millis: Given<'_>) -> PyResult<Self> {
                let inner = match &millis.0 {
                    None => Duration::ZERO,
                    Some(obj) => match obj.cast::<$dur>() {
                        Ok(same) => same.get().inner,
                        Err(_) => duration_from_py::<$n>(obj, "millis")?,
                    },
                };
                Ok(Self { inner })
            }

            /// Width of the magnitude in bits.
            #[classattr]
            const BITS: u32 = $n * 64;

            /// No time at all.
            #[classattr]
            #[pyo3(name = "ZERO")]
            fn zero() -> Self {
                Self {
                    inner: Duration::ZERO,
                }
            }

            /// Parse decimal digits with an optional leading minus sign.
            #[staticmethod]
            fn from_decimal(text: &str) -> PyResult<Self> {
                Ok(Self {
                    inner: Duration::from_decimal(text).py()?,
                })
            }
            /// The light-travel time across `distance_m` whole metres, rounded up to a whole tick, reading a
            /// tick as a millisecond. Exact. For another tick or grid use `Clock.light_ticks`.
            #[staticmethod]
            fn light_travel_time(distance_m: &Bound<'_, PyAny>) -> PyResult<Self> {
                let distance = uint_from_py::<$n>(distance_m, "distance_m")?;
                Ok(Self {
                    inner: Duration::light_travel_time_metres(distance).py()?,
                })
            }

            /// Whether the duration is below zero.
            #[getter]
            fn is_negative(&self) -> bool {
                self.inner.is_negative()
            }
            /// The duration in seconds, as a float, reading a tick as a millisecond. Raises
            /// `srkit.PrecisionLoss` beyond 2**53 ticks (285,000 years of milliseconds), where a float can
            /// no longer hold every tick. For ticks of another length use `Clock.local_seconds`.
            fn local_seconds(&self) -> PyResult<f64> {
                self.inner.local_seconds().py()
            }

            fn __int__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
                duration_to_py(py, self.inner)
            }
            fn __str__(&self) -> String {
                self.inner.to_string()
            }
            fn __repr__(&self) -> String {
                format!("{}({})", $dname, self.inner)
            }
            fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
                duration_to_py(py, self.inner)?.hash()
            }
            fn __bool__(&self) -> bool {
                self.inner != Duration::ZERO
            }
            fn __richcmp__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
                op: CompareOp,
            ) -> PyResult<Bound<'py, PyAny>> {
                if let Ok(d) = other.cast::<$dur>() {
                    let result = op.matches(compare(self.inner, d.get().inner));
                    Ok(result.into_pyobject(py)?.to_owned().into_any())
                } else if other.is_instance_of::<PyInt>() {
                    // An int may be wider than the duration, so compare as Python ints.
                    duration_to_py(py, self.inner)?.rich_compare(other, op)
                } else {
                    Ok(py.NotImplemented().into_bound(py))
                }
            }

            fn __neg__(&self) -> Self {
                Self {
                    inner: self.inner.negated(),
                }
            }
            fn __pos__(&self) -> Self {
                *self
            }
            fn __abs__(&self) -> Self {
                Self {
                    inner: Duration::from_parts(false, self.inner.magnitude_millis()),
                }
            }
            fn __add__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                match Self::arg(other)? {
                    Some(d) => {
                        let inner = self.inner.checked_add(d).py()?;
                        Ok(Bound::new(py, Self { inner })?.into_any())
                    }
                    None => Ok(py.NotImplemented().into_bound(py)),
                }
            }
            fn __radd__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                self.__add__(py, other)
            }
            fn __sub__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                match Self::arg(other)? {
                    Some(d) => {
                        let inner = self.inner.checked_sub(d).py()?;
                        Ok(Bound::new(py, Self { inner })?.into_any())
                    }
                    None => Ok(py.NotImplemented().into_bound(py)),
                }
            }
            fn __rsub__<'py>(
                &self,
                py: Python<'py>,
                other: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                match Self::arg(other)? {
                    Some(d) => {
                        let inner = d.checked_sub(self.inner).py()?;
                        Ok(Bound::new(py, Self { inner })?.into_any())
                    }
                    None => Ok(py.NotImplemented().into_bound(py)),
                }
            }
        }
    };
}

exact_time!(Time384, Duration384, 6, "Time384", "Duration384");
exact_time!(Time512, Duration512, 8, "Time512", "Duration512");
exact_time!(Time768, Duration768, 12, "Time768", "Duration768");
exact_time!(Time1024, Duration1024, 16, "Time1024", "Duration1024");
