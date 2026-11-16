//! Python bindings for SRKit: the extension module `srkit._srkit`.
//!
//! The Python package `srkit` (in `python/srkit/`) re-exports everything defined here and adds the vectorised
//! worldline methods. Each module below binds one part of `srkit-core`:
//!
//! | Module | Python names |
//! |---|---|
//! | [`segment`] | `Segment`, `TauState`, `RindlerHorizon` |
//! | [`worldline`] | `Worldline` (native base), `Radar` |
//! | [`trip`] | `Trip`, `Leg`, `TripState`, `FlipBurn`, `reach`, `flip_burn_tau`, `flip_burn_distance` |
//! | [`orbit`] | `CircularOrbit`, `OrbitState` |
//! | [`time`] | `Time384`, `Time512`, `Duration384`, `Duration512` |
//! | [`units`], [`plan`], [`causal`], [`functions`] | the submodules `srkit.units`, `srkit.plan`, `srkit.causal`, `srkit.kinematics`, `srkit.optics`, `srkit.rocket`, `srkit.kepler` |
//!
//! Numbers pass through unchanged: the bindings validate arguments but never normalise or round them, so a call
//! gives the same bits as the same call on the Rust core.
#![forbid(unsafe_code)]
// Comparisons are written `!(x > 0.0)` on purpose, as in the core: unlike `x <= 0.0` they also reject a NaN.
#![allow(clippy::neg_cmp_op_on_partial_ord, clippy::nonminimal_bool)]

use pyo3::prelude::*;

mod causal;
mod convert;
mod functions;
mod orbit;
mod plan;
mod segment;
mod time;
mod trip;
mod units;
mod worldline;

pyo3::create_exception!(
    srkit,
    NoIntersection,
    pyo3::exceptions::PyException,
    "A light cone does not meet a worldline where it was asked to.\n\n\
     Raised by a segment's `recv` and `send` when the light does not arrive within the range given, and by \
     `Worldline.seen_from` and `Worldline.light_cone_coordinates` when it never does. `Worldline.recv` and \
     `Worldline.send` do not raise: they return `inf` and `-inf`."
);

pyo3::create_exception!(
    srkit,
    PrecisionLoss,
    pyo3::exceptions::PyArithmeticError,
    "The answer cannot be given to the precision the function promises.\n\n\
     Raised when an exact difference of times is too large to be a float without losing milliseconds (beyond \
     2**53 ms), when a Doppler factor underflows, and when `causal.settle_arrival` is given an estimate \
     further off than its guard. A subclass of `ArithmeticError`, as `OverflowError` is."
);

/// The extension module behind the `srkit` package. Import `srkit` instead.
#[pymodule]
fn _srkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("NoIntersection", py.get_type::<NoIntersection>())?;
    m.add("PrecisionLoss", py.get_type::<PrecisionLoss>())?;

    m.add_class::<segment::Segment>()?;
    m.add_class::<segment::TauState>()?;
    m.add_class::<segment::RindlerHorizon>()?;
    m.add_class::<worldline::Worldline>()?;
    m.add_class::<worldline::Radar>()?;
    m.add_class::<trip::Trip>()?;
    m.add_class::<trip::Leg>()?;
    m.add_class::<trip::TripState>()?;
    m.add_class::<trip::FlipBurn>()?;
    m.add_class::<orbit::CircularOrbit>()?;
    m.add_class::<orbit::OrbitState>()?;
    m.add_class::<time::Time384>()?;
    m.add_class::<time::Time512>()?;
    m.add_class::<time::Duration384>()?;
    m.add_class::<time::Duration512>()?;
    m.add_function(wrap_pyfunction!(trip::reach, m)?)?;
    m.add_function(wrap_pyfunction!(trip::flip_burn_tau, m)?)?;
    m.add_function(wrap_pyfunction!(trip::flip_burn_distance, m)?)?;

    let kinematics = functions::kinematics_module(py)?;
    let optics = functions::optics_module(py)?;
    // `setattr`, not `add`: the names are already in the module's `__all__`.
    for name in functions::SHARED_WITH_OPTICS {
        optics.setattr(*name, kinematics.getattr(*name)?)?;
    }
    m.add("units", units::module(py)?)?;
    m.add("plan", plan::module(py)?)?;
    m.add("causal", causal::module(py)?)?;
    m.add("kinematics", kinematics)?;
    m.add("optics", optics)?;
    m.add("rocket", functions::rocket_module(py)?)?;
    m.add("kepler", functions::kepler_module(py)?)?;
    Ok(())
}
