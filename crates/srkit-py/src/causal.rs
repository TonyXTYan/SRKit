//! `srkit.causal`: exact causal decisions on stored integer times.
//!
//! Times are whole milliseconds (`Time384`, `Time512`) and displacements are whole metres as Python `int`s,
//! checked against the `i128` range. No float is involved in any decision.

use pyo3::{
    exceptions::{PyOverflowError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyString},
};
use srkit_core::{
    causal::{self as core, Causal, GridVec},
    time::Timestamp,
    Error,
};

use crate::{
    convert::{err, type_name, IntoPy},
    time::{duration_from_py, uint_to_py, Duration384, Duration512, Time384, Time512},
};

/// A displacement argument: a sequence of exactly three `int`s, in metres.
pub(crate) struct Grid(GridVec);

impl<'a, 'py> FromPyObject<'a, 'py> for Grid {
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if obj.is_instance_of::<PyString>() || obj.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "expected a sequence of three ints (metres), not text",
            ));
        }
        let len = obj.len().map_err(|_| {
            PyTypeError::new_err(format!(
                "expected a sequence of three ints (metres), got {}",
                type_name(&obj)
            ))
        })?;
        if len != 3 {
            return Err(PyValueError::new_err(format!(
                "expected three components (x, y, z) in metres, got {len}"
            )));
        }
        // The extraction goes through `__index__`: an int or a NumPy integer passes, a float does not.
        let component = |i: usize| -> PyResult<i128> {
            let item = obj.get_item(i)?;
            item.extract::<i128>().map_err(|e| {
                if e.is_instance_of::<PyTypeError>(obj.py()) {
                    PyTypeError::new_err(format!(
                        "a displacement is whole metres as ints, never floats; got {}",
                        type_name(&item)
                    ))
                } else {
                    PyOverflowError::new_err(
                        "a displacement component must fit a signed 128-bit integer of metres",
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

fn mixed_widths() -> PyErr {
    PyTypeError::new_err("the two times must have the same width (both Time384 or both Time512)")
}

fn not_a_time(obj: &Bound<'_, PyAny>) -> PyErr {
    PyTypeError::new_err(format!(
        "expected a Time384 or a Time512, got {}",
        type_name(obj)
    ))
}

/// The light-travel time over the displacement `dx` (three ints, metres), rounded up to whole milliseconds:
/// the smallest `n` with `c n >= 1000 |dx|`. Exact, by an integer square root. Returns an `int`.
///
/// Exactly one light-second (299,792,458 m) is 1000 ticks; one metre more is 1001.
#[pyfunction]
fn light_ticks(py: Python<'_>, dx: Grid) -> PyResult<Bound<'_, PyAny>> {
    let ticks = core::light_ticks::<8>(dx.0).py()?;
    uint_to_py(py, ticks.magnitude_millis())
}

/// The first tick at which a signal emitted at `emission` has reached a fixed point `dx` away (three ints,
/// metres). Returns a time of the same width as `emission`.
#[pyfunction]
fn first_arrival<'py>(
    py: Python<'py>,
    emission: &Bound<'py, PyAny>,
    dx: Grid,
) -> PyResult<Bound<'py, PyAny>> {
    if let Ok(t) = emission.cast::<Time384>() {
        let inner = core::first_arrival(t.get().inner, dx.0).py()?;
        Ok(Bound::new(py, Time384 { inner })?.into_any())
    } else if let Ok(t) = emission.cast::<Time512>() {
        let inner = core::first_arrival(t.get().inner, dx.0).py()?;
        Ok(Bound::new(py, Time512 { inner })?.into_any())
    } else {
        Err(not_a_time(emission))
    }
}

/// Whether a signal emitted at `emission` has reached a point `dx` away (three ints, metres) by the tick
/// `now`. Exact: true exactly from `first_arrival(emission, dx)` onwards, and false before the emission.
#[pyfunction]
fn arrived(emission: &Bound<'_, PyAny>, now: &Bound<'_, PyAny>, dx: Grid) -> PyResult<bool> {
    if let Ok(e) = emission.cast::<Time384>() {
        let n = now.cast::<Time384>().map_err(|_| mixed_widths())?;
        Ok(core::arrived(e.get().inner, n.get().inner, dx.0))
    } else if let Ok(e) = emission.cast::<Time512>() {
        let n = now.cast::<Time512>().map_err(|_| mixed_widths())?;
        Ok(core::arrived(e.get().inner, n.get().inner, dx.0))
    } else {
        Err(not_a_time(emission))
    }
}

/// Classify the separation between two events `dt` apart in time (an `int` of milliseconds or a duration; the
/// sign does not matter) and `dx` apart in space (three ints, metres).
///
/// Returns 'timelike' (the later event is strictly inside the earlier one's light cone), 'null' (light from
/// the earlier reaches the later exactly) or 'spacelike' (neither can influence the other). Exact.
#[pyfunction]
fn classify(dt: &Bound<'_, PyAny>, dx: Grid) -> PyResult<&'static str> {
    let causal = if let Ok(d) = dt.cast::<Duration384>() {
        core::classify(d.get().inner, dx.0)
    } else if let Ok(d) = dt.cast::<Duration512>() {
        core::classify(d.get().inner, dx.0)
    } else {
        core::classify(duration_from_py::<8>(dt, "dt")?, dx.0)
    };
    Ok(match causal {
        Causal::Timelike => "timelike",
        Causal::Null => "null",
        Causal::Spacelike => "spacelike",
    })
}

/// Round a displacement in metres, given as floats, to the one-metre grid (to nearest, ties away from zero).
/// Returns three ints. Raises `OverflowError` for a component that is not finite or is beyond the grid.
#[pyfunction]
fn grid_from_metres(x: f64, y: f64, z: f64) -> PyResult<(i128, i128, i128)> {
    let g = GridVec::from_metres(x, y, z).py()?;
    Ok((g.x, g.y, g.z))
}

/// Run the core's search with a Python callback, keeping the callback's own exception if it raises.
fn settle<const N: usize>(
    emission: Timestamp<N>,
    estimate: Timestamp<N>,
    guard: u32,
    mut displacement: impl FnMut(Timestamp<N>) -> PyResult<GridVec>,
) -> PyResult<Timestamp<N>> {
    let mut raised: Option<PyErr> = None;
    let result = core::settle_arrival(emission, estimate, guard, |tick| {
        displacement(tick).map_err(|e| {
            raised = Some(e);
            Error::InvalidInput
        })
    });
    match (result, raised) {
        (_, Some(e)) => Err(e),
        (Ok(tick), None) => Ok(tick),
        (Err(e), None) => Err(err(e)),
    }
}

/// The first tick at which a signal reaches a moving receiver, settled exactly from a floating-point estimate.
///
/// `displacement_at(tick)` is called with a time and must return the receiver's position relative to the
/// emission point at that tick, as three ints in metres (evaluate the receiver's worldline in floats and round
/// with `grid_from_metres`). Starting from `estimate`, the search steps back while the signal has already
/// arrived and forward until it has, so the result is the first tick that passes the exact check and can never
/// be a tick early. `guard` bounds the number of steps each way; an estimate further off than that raises
/// `srkit.PrecisionLoss`. Returns a time of the same width as `emission`.
#[pyfunction]
#[pyo3(signature = (emission, estimate, displacement_at, guard=8))]
fn settle_arrival<'py>(
    py: Python<'py>,
    emission: &Bound<'py, PyAny>,
    estimate: &Bound<'py, PyAny>,
    displacement_at: &Bound<'py, PyAny>,
    guard: u32,
) -> PyResult<Bound<'py, PyAny>> {
    if let Ok(e) = emission.cast::<Time384>() {
        let est = estimate.cast::<Time384>().map_err(|_| mixed_widths())?;
        let inner = settle(e.get().inner, est.get().inner, guard, |tick| {
            let tick = Bound::new(py, Time384 { inner: tick })?;
            Ok(displacement_at.call1((tick,))?.extract::<Grid>()?.0)
        })?;
        Ok(Bound::new(py, Time384 { inner })?.into_any())
    } else if let Ok(e) = emission.cast::<Time512>() {
        let est = estimate.cast::<Time512>().map_err(|_| mixed_widths())?;
        let inner = settle(e.get().inner, est.get().inner, guard, |tick| {
            let tick = Bound::new(py, Time512 { inner: tick })?;
            Ok(displacement_at.call1((tick,))?.extract::<Grid>()?.0)
        })?;
        Ok(Bound::new(py, Time512 { inner })?.into_any())
    } else {
        Err(not_a_time(emission))
    }
}

const DOC: &str = "\
Exact causal decisions on stored integer times.

Geometry in srkit is floating point, which is right for drawing and wrong for deciding. \"Has this signal
arrived by tick n?\" and \"is event B inside A's light cone?\" are answered here in integer arithmetic: times
are whole milliseconds (`srkit.Time384`, `srkit.Time512`) and displacements are whole metres, given as three
Python ints (each within the range of a signed 128-bit integer, 1.8e22 light-years). The speed of light is
exactly 299,792,458 m/s, so

    arrived  <=>  c * dt_ms >= 1000 * |dx_m|

holds or fails with no rounding at all. Compute where things are in floating point, round the endpoints to the
grid (`grid_from_metres`), then decide exactly. Floats are refused as displacements.
";

pub(crate) fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.causal")?;
    m.setattr("__doc__", DOC)?;
    m.add("C_METRES_PER_SECOND", core::C_METRES_PER_SECOND)?;
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
