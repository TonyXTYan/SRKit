//! Shared conversions: three-vectors, error mapping, argument checks and `repr` text.

use pyo3::{
    exceptions::{PyOverflowError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyString, PyTuple},
};
use srkit_core::{kinematics::Event, Error, Vec3};

use crate::{NoIntersection, PrecisionLoss};

/// A Python three-tuple of floats.
pub(crate) type Tup3 = (f64, f64, f64);

/// A three-vector argument: any sequence of exactly three real numbers (a tuple, a list, a NumPy array).
#[derive(Clone, Copy, Debug)]
pub(crate) struct V3(pub Vec3);

impl<'a, 'py> FromPyObject<'a, 'py> for V3 {
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(t) = obj.cast::<PyTuple>() {
            if t.len() == 3 {
                if let Ok((x, y, z)) = t.extract::<Tup3>() {
                    return Ok(Self(Vec3::new(x, y, z)));
                }
            }
        }
        if obj.is_instance_of::<PyString>() || obj.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "expected a sequence of three numbers, not text",
            ));
        }
        let len = obj.len().map_err(|_| {
            PyTypeError::new_err(format!(
                "expected a sequence of three numbers, got {}",
                type_name(&obj)
            ))
        })?;
        if len != 3 {
            return Err(PyValueError::new_err(format!(
                "expected three components (x, y, z), got {len}"
            )));
        }
        let component = |i: usize| -> PyResult<f64> { obj.get_item(i)?.extract::<f64>() };
        Ok(Self(Vec3::new(component(0)?, component(1)?, component(2)?)))
    }
}

/// The name of an object's type, for error messages.
pub(crate) fn type_name(obj: &Bound<'_, PyAny>) -> String {
    obj.get_type()
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|_| "object".to_owned())
}

/// A vector as a Python tuple.
pub(crate) fn tup(v: Vec3) -> Tup3 {
    (v.x, v.y, v.z)
}

/// A vector from a stored Python tuple.
pub(crate) fn vec(t: Tup3) -> Vec3 {
    Vec3::new(t.0, t.1, t.2)
}

/// An event from its time and position arguments.
pub(crate) fn event(t: f64, x: V3) -> Event {
    Event::new(t, x.0)
}

/// Map a core error to the Python exception that stands for it.
///
/// `InvalidInput`, `OutOfRange` and `Discontinuous` are `ValueError`; `Overflow` is `OverflowError`;
/// `NoIntersection` and `PrecisionLoss` have exceptions of their own.
pub(crate) fn err(e: Error) -> PyErr {
    let message = e.to_string();
    match e {
        Error::InvalidInput | Error::OutOfRange | Error::Discontinuous => {
            PyValueError::new_err(message)
        }
        Error::Overflow => PyOverflowError::new_err(message),
        Error::NoIntersection => NoIntersection::new_err(message),
        Error::PrecisionLoss => PrecisionLoss::new_err(message),
    }
}

/// `Result<T, srkit_core::Error>` to `PyResult<T>`.
pub(crate) trait IntoPy<T> {
    fn py(self) -> PyResult<T>;
}

impl<T> IntoPy<T> for Result<T, Error> {
    fn py(self) -> PyResult<T> {
        self.map_err(err)
    }
}

/// A `ValueError` with the given message.
pub(crate) fn invalid(message: impl Into<String>) -> PyErr {
    PyValueError::new_err(message.into())
}

/// Reject a NaN or an infinity.
pub(crate) fn finite(name: &str, x: f64) -> PyResult<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid(format!("{name} must be finite, got {}", fl(x))))
    }
}

/// Reject a value that is not a finite number above zero.
pub(crate) fn positive(name: &str, x: f64) -> PyResult<f64> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(invalid(format!(
            "{name} must be finite and positive, got {}",
            fl(x)
        )))
    }
}

/// Reject a vector with a NaN or an infinite component.
pub(crate) fn finite_vec(name: &str, v: Vec3) -> PyResult<Vec3> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid(format!("{name} must have finite components")))
    }
}

/// How far a direction's length may be from 1 and still count as a unit vector.
const UNIT_TOLERANCE: f64 = 1e-9;

/// Reject a direction that is not a unit vector. The vector is passed through unchanged, so that the same
/// numbers give the same bits as in the Rust core.
pub(crate) fn unit_vec(name: &str, v: Vec3) -> PyResult<Vec3> {
    finite_vec(name, v)?;
    let norm = v.norm();
    if (norm - 1.0).abs() <= UNIT_TOLERANCE {
        Ok(v)
    } else {
        Err(invalid(format!(
            "{name} must be a unit vector, its length is {}",
            fl(norm)
        )))
    }
}

/// Reject a velocity at or above the speed of light.
pub(crate) fn subluminal(name: &str, v: Vec3) -> PyResult<Vec3> {
    finite_vec(name, v)?;
    let speed = v.norm();
    if speed < 1.0 {
        Ok(v)
    } else {
        Err(invalid(format!(
            "{name} must be slower than light (a fraction of c below 1), its speed is {}",
            fl(speed)
        )))
    }
}

/// A float as Python writes it, near enough: `1.0`, `1e100`, `inf`, `nan`.
pub(crate) fn fl(x: f64) -> String {
    if x.is_nan() {
        "nan".to_owned()
    } else if x.is_infinite() {
        if x > 0.0 { "inf" } else { "-inf" }.to_owned()
    } else {
        format!("{x:?}")
    }
}

/// The text a value has inside a `repr`.
pub(crate) trait Repr {
    fn repr(&self, py: Python<'_>) -> String;
}

impl Repr for f64 {
    fn repr(&self, _: Python<'_>) -> String {
        fl(*self)
    }
}
impl Repr for bool {
    fn repr(&self, _: Python<'_>) -> String {
        if *self { "True" } else { "False" }.to_owned()
    }
}
impl Repr for &'static str {
    fn repr(&self, _: Python<'_>) -> String {
        format!("'{self}'")
    }
}
impl Repr for Tup3 {
    fn repr(&self, _: Python<'_>) -> String {
        format!("({}, {}, {})", fl(self.0), fl(self.1), fl(self.2))
    }
}
impl Repr for Vec3 {
    fn repr(&self, py: Python<'_>) -> String {
        tup(*self).repr(py)
    }
}
impl<T: Repr> Repr for Option<T> {
    fn repr(&self, py: Python<'_>) -> String {
        match self {
            Some(value) => value.repr(py),
            None => "None".to_owned(),
        }
    }
}
impl<T> Repr for Py<T> {
    fn repr(&self, py: Python<'_>) -> String {
        self.bind(py)
            .as_any()
            .repr()
            .map(|r| r.to_string())
            .unwrap_or_else(|_| "<?>".to_owned())
    }
}

/// `Name(field=value, ...)` from the named fields of `self`.
macro_rules! repr_fields {
    ($py:expr, $name:expr, $this:ident, $($field:ident),+ $(,)?) => {{
        let parts = [$(format!(
            "{}={}",
            stringify!($field),
            $crate::convert::Repr::repr(&$this.$field, $py)
        )),+];
        format!("{}({})", $name, parts.join(", "))
    }};
}
pub(crate) use repr_fields;

/// `Name(field=value, ...)` from explicit values.
macro_rules! repr_values {
    ($py:expr, $name:expr, $($field:ident = $value:expr),+ $(,)?) => {{
        let parts = [$(format!(
            "{}={}",
            stringify!($field),
            $crate::convert::Repr::repr(&$value, $py)
        )),+];
        format!("{}({})", $name, parts.join(", "))
    }};
}
pub(crate) use repr_values;
