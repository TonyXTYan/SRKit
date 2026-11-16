//! `srkit.plan`: itineraries, as sequences of segments that join continuously.

use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple};
use srkit_core::{plan as core, Vec3};

use crate::{
    convert::{
        finite, finite_vec, fl, invalid, positive, repr_fields, repr_values, subluminal, tup,
        type_name, unit_vec, vec, Tup3, V3,
    },
    orbit::{orbit_of_speed, CircularOrbit},
    segment::Segment,
};

fn phase_name(phase: core::Phase) -> &'static str {
    match phase {
        core::Phase::Rest => "rest",
        core::Phase::Boost => "boost",
        core::Phase::Coast => "coast",
        core::Phase::Brake => "brake",
        core::Phase::Orbit => "orbit",
        core::Phase::Jump => "jump",
    }
}

fn phase_of(name: &str) -> PyResult<core::Phase> {
    Ok(match name {
        "rest" => core::Phase::Rest,
        "boost" => core::Phase::Boost,
        "coast" => core::Phase::Coast,
        "brake" => core::Phase::Brake,
        "orbit" => core::Phase::Orbit,
        "jump" => core::Phase::Jump,
        other => {
            return Err(invalid(format!(
                "phase must be 'rest', 'boost', 'coast', 'brake', 'orbit' or 'jump', got '{other}'"
            )))
        }
    })
}

/// Where and when a plan (or its next piece) begins.
///
/// `Start(t=0.0, tau=0.0, x=(0.0, 0.0, 0.0))`: the coordinate time, the traveller's proper time and the
/// position.
#[pyclass(frozen, get_all, module = "srkit.plan", name = "Start")]
pub struct Start {
    /// Coordinate time.
    t: f64,
    /// The traveller's proper time.
    tau: f64,
    /// Position `(x, y, z)`.
    x: Tup3,
}

#[pymethods]
impl Start {
    #[new]
    #[pyo3(signature = (t=0.0, tau=0.0, x=None))]
    fn new(t: f64, tau: f64, x: Option<V3>) -> PyResult<Self> {
        Ok(Self {
            t: finite("t", t)?,
            tau: finite("tau", tau)?,
            x: tup(finite_vec("x", x.map_or(Vec3::ZERO, |v| v.0))?),
        })
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "Start", self, t, tau, x)
    }
}

/// The state at the end of a piece or a plan.
#[pyclass(frozen, get_all, module = "srkit.plan", name = "End")]
pub struct End {
    /// Coordinate time.
    t: f64,
    /// The traveller's proper time.
    tau: f64,
    /// Position `(x, y, z)`.
    x: Tup3,
    /// Velocity as a fraction of c.
    v: Tup3,
    /// Rapidity along the line of the last piece.
    phi: f64,
    /// Proper acceleration at the end of the last piece.
    acc: f64,
}

impl From<core::End> for End {
    fn from(e: core::End) -> Self {
        Self {
            t: e.t,
            tau: e.tau,
            x: tup(e.x),
            v: tup(e.v),
            phi: e.phi,
            acc: e.acc,
        }
    }
}

#[pymethods]
impl End {
    /// A `Start` at this time, proper time and position, to carry on from here.
    fn start(&self) -> Start {
        Start {
            t: self.t,
            tau: self.tau,
            x: self.x,
        }
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "End", self, t, tau, x, v, phi, acc)
    }
}

/// A start argument: a `Start`, or the `End` of an earlier piece.
pub(crate) struct StartArg(core::Start);

impl<'a, 'py> FromPyObject<'a, 'py> for StartArg {
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        let (t, tau, x) = if let Ok(s) = obj.cast::<Start>() {
            let s = s.get();
            (s.t, s.tau, s.x)
        } else if let Ok(e) = obj.cast::<End>() {
            let e = e.get();
            (e.t, e.tau, e.x)
        } else {
            return Err(PyTypeError::new_err(format!(
                "expected a plan.Start (or the plan.End of an earlier piece), got {}",
                type_name(&obj)
            )));
        };
        Ok(Self(core::Start { t, tau, x: vec(x) }))
    }
}

/// One piece of a plan: a segment, the proper time at its start, how much proper time it takes and what it is
/// for.
#[pyclass(frozen, module = "srkit.plan", name = "Piece")]
pub struct Piece {
    inner: core::Piece,
}

impl From<core::Piece> for Piece {
    fn from(inner: core::Piece) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Piece {
    /// The piece's motion, as a `Segment`.
    #[getter]
    fn segment(&self) -> Segment {
        self.inner.seg.into()
    }
    /// The traveller's proper time at the start of the piece.
    #[getter]
    fn tau0(&self) -> f64 {
        self.inner.tau0
    }
    /// Proper time the piece takes.
    #[getter]
    fn dtau(&self) -> f64 {
        self.inner.dtau
    }
    /// What the piece is for: 'rest', 'boost', 'coast', 'brake', 'orbit' or 'jump'.
    #[getter]
    fn phase(&self) -> &'static str {
        phase_name(self.inner.phase)
    }
    /// The state at the end of the piece.
    fn end(&self) -> End {
        self.inner.end().into()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        let segment = Py::new(py, self.segment()).ok();
        repr_values!(
            py,
            "Piece",
            segment = segment,
            tau0 = self.inner.tau0,
            dtau = self.inner.dtau,
            phase = phase_name(self.inner.phase)
        )
    }
}

fn piece_tuple<'py>(py: Python<'py>, pieces: &core::Pieces) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        pieces.iter().map(|p| Piece::from(*p)).collect::<Vec<_>>(),
    )
}

fn segment_tuple<'py>(py: Python<'py>, pieces: &core::Pieces) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        pieces
            .iter()
            .map(|p| Segment::from(p.seg))
            .collect::<Vec<_>>(),
    )
}

/// A straight run and its summary: what `line` returns.
#[pyclass(frozen, module = "srkit.plan", name = "Line")]
pub struct Line {
    inner: core::Line,
}

#[pymethods]
impl Line {
    /// The pieces in order, as a tuple of `Piece`: a boost, a coast if the cap is reached, a brake. Empty for
    /// a run of no length.
    #[getter]
    fn pieces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        piece_tuple(py, &self.inner.pieces)
    }
    /// The pieces' segments in order, as a tuple of `Segment`, ready for `Worldline`.
    #[getter]
    fn segments<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        segment_tuple(py, &self.inner.pieces)
    }
    /// The state at the end of the run. Its position is pinned exactly on the target.
    #[getter]
    fn end(&self) -> End {
        self.inner.end.into()
    }
    /// Unit vector along the run.
    #[getter]
    fn n(&self) -> Tup3 {
        tup(self.inner.n)
    }
    /// Total change of rapidity the drive supplies (boost plus brake): what the rocket equation charges for.
    #[getter]
    fn dphi(&self) -> f64 {
        self.inner.dphi
    }
    /// Length of the run.
    #[getter]
    fn length(&self) -> f64 {
        self.inner.length
    }
    /// Peak rapidity.
    #[getter]
    fn phi_peak(&self) -> f64 {
        self.inner.phi_peak
    }
    /// Length coasted at the cap (zero if the cap is not reached).
    #[getter]
    fn coast_length(&self) -> f64 {
        self.inner.coast_length
    }
    fn __repr__(&self) -> String {
        format!(
            "Line({} piece{}, length={}, dphi={}, phi_peak={}, coast_length={}, end.t={}, end.tau={})",
            self.inner.pieces.len(),
            if self.inner.pieces.len() == 1 { "" } else { "s" },
            fl(self.inner.length),
            fl(self.inner.dphi),
            fl(self.inner.phi_peak),
            fl(self.inner.coast_length),
            fl(self.inner.end.t),
            fl(self.inner.end.tau)
        )
    }
}

/// A braking burn to rest and its summary: what `brake` and `orbit_exit` return.
#[pyclass(frozen, module = "srkit.plan", name = "Brake")]
pub struct Brake {
    piece: Option<core::Piece>,
    end: core::End,
    dphi: f64,
}

impl From<(Option<core::Piece>, core::End, f64)> for Brake {
    fn from((piece, end, dphi): (Option<core::Piece>, core::End, f64)) -> Self {
        Self { piece, end, dphi }
    }
}

#[pymethods]
impl Brake {
    /// The braking piece in a tuple, or an empty tuple if the ship was already at rest.
    #[getter]
    fn pieces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.piece.map(Piece::from))
    }
    /// The braking segment in a tuple, or an empty tuple if the ship was already at rest.
    #[getter]
    fn segments<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.piece.map(|p| Segment::from(p.seg)))
    }
    /// The state at rest.
    #[getter]
    fn end(&self) -> End {
        self.end.into()
    }
    /// The rapidity spent.
    #[getter]
    fn dphi(&self) -> f64 {
        self.dphi
    }
    fn __repr__(&self) -> String {
        format!(
            "Brake({} piece{}, dphi={}, end.t={}, end.tau={})",
            usize::from(self.piece.is_some()),
            if self.piece.is_some() { "" } else { "s" },
            fl(self.dphi),
            fl(self.end.t),
            fl(self.end.tau)
        )
    }
}

/// An orthonormal pair `(e1, e2)` and their normal: what `plane_through` returns.
#[pyclass(frozen, get_all, module = "srkit.plan", name = "Plane")]
pub struct Plane {
    /// Unit vector from the centre to the point.
    e1: Tup3,
    /// Unit vector at right angles to `e1` in the plane.
    e2: Tup3,
    /// Unit normal of the plane.
    normal: Tup3,
}

#[pymethods]
impl Plane {
    fn __repr__(&self, py: Python<'_>) -> String {
        repr_fields!(py, "Plane", self, e1, e2, normal)
    }
}

/// The elements of a circle to fly: centre, radius, angular velocity, plane and the angle at the join.
///
/// `OrbitElements(centre, r, omega, e1=(1, 0, 0), e2=(0, 1, 0), a0=0.0)`. `orbit_entry` returns one for the
/// circle it joins.
#[pyclass(frozen, module = "srkit.plan", name = "OrbitElements")]
pub struct OrbitElements {
    inner: core::OrbitElements,
}

#[pymethods]
impl OrbitElements {
    #[new]
    #[pyo3(signature = (centre, r, omega, e1=None, e2=None, a0=0.0))]
    fn new(
        centre: V3,
        r: f64,
        omega: f64,
        e1: Option<V3>,
        e2: Option<V3>,
        a0: f64,
    ) -> PyResult<Self> {
        let e1 = unit_vec("e1", e1.map_or(Vec3::X, |v| v.0))?;
        let e2 = unit_vec("e2", e2.map_or(Vec3::Y, |v| v.0))?;
        if e1.dot(e2).abs() > 1e-9 {
            return Err(invalid("e1 and e2 must be at right angles"));
        }
        finite("omega", omega)?;
        let orbit = orbit_of_speed(r, (r * omega).abs())?;
        Ok(Self {
            inner: core::OrbitElements {
                centre: finite_vec("centre", centre.0)?,
                r,
                omega,
                e1,
                e2,
                a0: finite("a0", a0)?,
                orbit,
            },
        })
    }
    /// Centre of the circle.
    #[getter]
    fn centre(&self) -> Tup3 {
        tup(self.inner.centre)
    }
    /// Radius.
    #[getter]
    fn r(&self) -> f64 {
        self.inner.r
    }
    /// Angular velocity in coordinate time, radians per unit of time. Positive runs from `e1` towards `e2`.
    #[getter]
    fn omega(&self) -> f64 {
        self.inner.omega
    }
    /// First unit vector of the plane.
    #[getter]
    fn e1(&self) -> Tup3 {
        tup(self.inner.e1)
    }
    /// Second unit vector of the plane.
    #[getter]
    fn e2(&self) -> Tup3 {
        tup(self.inner.e2)
    }
    /// Angle from `e1` at the join, in radians.
    #[getter]
    fn a0(&self) -> f64 {
        self.inner.a0
    }
    /// The orbit's speed, period and clock rates, as a `CircularOrbit`.
    #[getter]
    fn orbit(&self) -> CircularOrbit {
        self.inner.orbit.into()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        let g = &self.inner;
        repr_values!(
            py,
            "OrbitElements",
            centre = g.centre,
            r = g.r,
            omega = g.omega,
            e1 = g.e1,
            e2 = g.e2,
            a0 = g.a0
        )
    }
}

/// The pieces that bring a ship onto a circular orbit, and the circle's elements at the join: what
/// `orbit_entry` returns.
#[pyclass(frozen, module = "srkit.plan", name = "OrbitEntry")]
pub struct OrbitEntry {
    inner: core::OrbitEntry,
}

#[pymethods]
impl OrbitEntry {
    /// The pieces of the approach in order, as a tuple of `Piece`. The circle itself is not among them: add it
    /// with `circ_piece`, so that it can be open-ended.
    #[getter]
    fn pieces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        piece_tuple(py, &self.inner.pieces)
    }
    /// The approach's segments in order, as a tuple of `Segment`.
    #[getter]
    fn segments<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        segment_tuple(py, &self.inner.pieces)
    }
    /// The state at the join, moving at the orbit's speed along its tangent.
    #[getter]
    fn end(&self) -> End {
        self.inner.end.into()
    }
    /// The circle's elements at the join.
    #[getter]
    fn orbit(&self) -> OrbitElements {
        OrbitElements {
            inner: self.inner.orb,
        }
    }
    /// Total change of rapidity the drive supplies.
    #[getter]
    fn dphi(&self) -> f64 {
        self.inner.dphi
    }
    /// True when the entry started near the centre (a run to a point behind the join, then a boost along the
    /// tangent); false for the tangent run from outside.
    #[getter(from_inside)]
    fn is_from_inside(&self) -> bool {
        self.inner.from_inside
    }
    fn __repr__(&self) -> String {
        format!(
            "OrbitEntry({} pieces, dphi={}, from_inside={}, end.t={}, end.tau={})",
            self.inner.pieces.len(),
            fl(self.inner.dphi),
            if self.inner.from_inside {
                "True"
            } else {
                "False"
            },
            fl(self.inner.end.t),
            fl(self.inner.end.tau)
        )
    }
}

/// A straight run from `start` to the point `x1`: boost at proper acceleration `a`, coast if the rapidity
/// reaches `phi_cap`, brake at `a_brake` (by default `a`).
///
/// The rapidity along the line is `phi0` at the start and `phi1` at the end (both zero by default: rest to
/// rest). `phi_cap` is a rapidity; to cap the Lorentz factor at `g` pass `math.acosh(g)`. `n_default` is the
/// direction to use when the run has no length but a change of speed is asked for.
///
/// Returns a `Line`. Raises `ValueError` if the run is too short to change speed as asked.
#[pyfunction]
#[pyo3(
    signature = (start, x1, a, phi_cap=f64::INFINITY, phi0=0.0, phi1=0.0, a_brake=None, n_default=None),
    text_signature = "(start, x1, a, phi_cap=1e999, phi0=0.0, phi1=0.0, a_brake=None, n_default=None)"
)]
#[allow(clippy::too_many_arguments)]
fn line(
    start: StartArg,
    x1: V3,
    a: f64,
    phi_cap: f64,
    phi0: f64,
    phi1: f64,
    a_brake: Option<f64>,
    n_default: Option<V3>,
) -> PyResult<Line> {
    positive("a", a)?;
    let a_brake = positive("a_brake", a_brake.unwrap_or(a))?;
    if !(phi_cap > 0.0) {
        return Err(invalid(format!(
            "phi_cap must be positive (it may be infinite), got {}",
            fl(phi_cap)
        )));
    }
    for (name, phi) in [("phi0", phi0), ("phi1", phi1)] {
        if !(phi.is_finite() && phi >= 0.0) {
            return Err(invalid(format!(
                "{name} must be finite and not negative, got {}",
                fl(phi)
            )));
        }
    }
    let n_default = n_default.map(|n| unit_vec("n_default", n.0)).transpose()?;
    core::line(
        start.0,
        finite_vec("x1", x1.0)?,
        a,
        phi_cap,
        phi0,
        phi1,
        a_brake,
        n_default,
    )
    .map(|inner| Line { inner })
    .ok_or_else(|| invalid("the run is too short for the change of speed asked for"))
}

/// Brake from rapidity `phi0` along the unit vector `n` to rest, at proper acceleration `a`.
///
/// Returns a `Brake`: the piece (none if already at rest), the end state and the rapidity spent.
#[pyfunction]
fn brake(start: StartArg, n: V3, phi0: f64, a: f64) -> PyResult<Brake> {
    positive("a", a)?;
    finite("phi0", phi0)?;
    Ok(core::brake(start.0, unit_vec("n", n.0)?, phi0, a).into())
}

/// Leave a circle while moving at velocity `v` (a fraction of c): brake along the tangent to rest at proper
/// acceleration `a`. Returns a `Brake`.
#[pyfunction]
fn orbit_exit(start: StartArg, v: V3, a: f64) -> PyResult<Brake> {
    positive("a", a)?;
    let v = subluminal("v", v.0)?;
    if v.norm() == 0.0 {
        return Ok(core::brake(start.0, Vec3::X, 0.0, a).into());
    }
    Ok(core::orbit_exit(start.0, v, a).into())
}

/// The plane through `centre` that holds the point `x` and is tilted least from the xy plane. `e1` points from
/// the centre to `x`; when `x` is at the centre, `e1_default` (by default the x axis) is used instead.
#[pyfunction]
#[pyo3(signature = (centre, x, e1_default=None))]
fn plane_through(centre: V3, x: V3, e1_default: Option<V3>) -> PyResult<Plane> {
    let e1_default = e1_default
        .map(|v| finite_vec("e1_default", v.0))
        .transpose()?;
    let p = core::plane_through(
        finite_vec("centre", centre.0)?,
        finite_vec("x", x.0)?,
        e1_default,
    );
    Ok(Plane {
        e1: tup(p.e1),
        e2: tup(p.e2),
        normal: tup(p.normal),
    })
}

/// Enter a circular orbit of radius `r` at speed `beta` about `centre`, from rest at `start`, at proper
/// acceleration `a`. The orbit runs counter-clockwise in `(e1, e2)`.
///
/// From outside, the run is the straight line tangent to the circle at its join (boost, coast, brake to the
/// orbit's speed), so the velocity is continuous. From near the centre: a rest-to-rest run to a point behind
/// the join, then a boost along the tangent. `phi_cap` caps the rapidity of the run. `e1` sets the plane when
/// the start is at the centre, `a0` the join angle when entering from near the centre, and `a_tangent` the
/// proper acceleration of the boost along the tangent when it differs from `a`.
///
/// Returns an `OrbitEntry`. Raises `ValueError` if a run is too short for its change of speed.
#[pyfunction]
#[pyo3(
    signature = (start, centre, r, beta, a, phi_cap=f64::INFINITY, *, e1=None, a0=None, a_tangent=None),
    text_signature = "(start, centre, r, beta, a, phi_cap=1e999, *, e1=None, a0=None, a_tangent=None)"
)]
#[allow(clippy::too_many_arguments)]
fn orbit_entry(
    start: StartArg,
    centre: V3,
    r: f64,
    beta: f64,
    a: f64,
    phi_cap: f64,
    e1: Option<V3>,
    a0: Option<f64>,
    a_tangent: Option<f64>,
) -> PyResult<OrbitEntry> {
    orbit_of_speed(r, beta)?;
    positive("a", a)?;
    if !(phi_cap > 0.0) {
        return Err(invalid(format!(
            "phi_cap must be positive (it may be infinite), got {}",
            fl(phi_cap)
        )));
    }
    let opts = core::EntryOptions {
        e1: e1.map(|v| unit_vec("e1", v.0)).transpose()?,
        a0: a0.map(|v| finite("a0", v)).transpose()?,
        a_tangent: a_tangent.map(|v| positive("a_tangent", v)).transpose()?,
    };
    core::orbit_entry(
        start.0,
        finite_vec("centre", centre.0)?,
        r,
        beta,
        a,
        phi_cap,
        opts,
    )
    .map(|inner| OrbitEntry { inner })
    .ok_or_else(|| invalid("a run is too short for its change of speed"))
}

/// Stay at rest for proper time `dtau` (which may be infinite).
#[pyfunction]
fn rest(start: StartArg, dtau: f64) -> PyResult<Piece> {
    if !(dtau >= 0.0) {
        return Err(invalid(format!(
            "dtau must not be negative, got {}",
            fl(dtau)
        )));
    }
    Ok(core::rest(start.0, dtau).into())
}

/// Accelerate at proper acceleration `alpha` (signed) along the unit vector `n`, from rapidity `phi0` to
/// `phi1` along `n`. `phase` labels the piece: 'boost' by default.
#[pyfunction]
#[pyo3(signature = (start, n, alpha, phi0, phi1, phase="boost"))]
fn hyp_piece(
    start: StartArg,
    n: V3,
    alpha: f64,
    phi0: f64,
    phi1: f64,
    phase: &str,
) -> PyResult<Piece> {
    if !(alpha.is_finite() && alpha != 0.0) {
        return Err(invalid(format!(
            "alpha must be finite and not zero, got {}",
            fl(alpha)
        )));
    }
    finite("phi0", phi0)?;
    finite("phi1", phi1)?;
    if (phi1 - phi0) / alpha < 0.0 {
        return Err(invalid(
            "alpha has the wrong sign to take the rapidity from phi0 to phi1",
        ));
    }
    Ok(core::hyp_piece(
        start.0,
        unit_vec("n", n.0)?,
        alpha,
        phi0,
        phi1,
        phase_of(phase)?,
    )
    .into())
}

/// Coast at velocity `v` (a fraction of c) for proper time `dtau`.
#[pyfunction]
fn lin_piece(start: StartArg, v: V3, dtau: f64) -> PyResult<Piece> {
    if !(dtau.is_finite() && dtau >= 0.0) {
        return Err(invalid(format!(
            "dtau must be finite and not negative, got {}",
            fl(dtau)
        )));
    }
    Ok(core::lin_piece(start.0, subluminal("v", v.0)?, dtau).into())
}

/// A light-like jump to the point `x1`: it takes the light-travel time in coordinate time and no proper time.
#[pyfunction]
fn null_piece(start: StartArg, x1: V3) -> PyResult<Piece> {
    Ok(core::null_piece(start.0, finite_vec("x1", x1.0)?).into())
}

/// Fly the circle `orbit` (an `OrbitElements`) for proper time `dtau`, by default for ever, starting at its
/// join angle.
#[pyfunction]
#[pyo3(
    signature = (start, orbit, dtau=f64::INFINITY),
    text_signature = "(start, orbit, dtau=1e999)"
)]
fn circ_piece(start: StartArg, orbit: &OrbitElements, dtau: f64) -> PyResult<Piece> {
    if !(dtau >= 0.0) {
        return Err(invalid(format!(
            "dtau must not be negative, got {}",
            fl(dtau)
        )));
    }
    Ok(core::circ_piece(start.0, &orbit.inner, dtau).into())
}

const DOC: &str = "\
Itineraries: sequences of segments that join continuously, with no instant change of velocity.

A plan starts from a `Start` (a time, a proper time and a position) and is built from pieces, each a
closed-form `Segment` with its proper duration. `line` is the workhorse: a straight run that boosts, coasts if
the speed is capped, and brakes, between any two rapidities along the line. `orbit_entry` joins a circular
orbit along a tangent so that the velocity is continuous at the join.

Results carry their pieces as `Piece` objects (`.pieces`) and as plain `Segment` objects (`.segments`), and an
`End` (`.end`) that the next call accepts as its start:

    out = plan.line(plan.Start(), (4.2465, 0, 0), a=units.G0_PER_YEAR)
    ship = srkit.Worldline([*out.segments, plan.rest(out.end, math.inf).segment])

Natural units, c = 1: accelerations are rates of change of rapidity per unit of time.
";

pub(crate) fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let m = PyModule::new(py, "srkit.plan")?;
    m.setattr("__doc__", DOC)?;
    m.add_class::<Start>()?;
    m.add_class::<End>()?;
    m.add_class::<Piece>()?;
    m.add_class::<Line>()?;
    m.add_class::<Brake>()?;
    m.add_class::<Plane>()?;
    m.add_class::<OrbitElements>()?;
    m.add_class::<OrbitEntry>()?;
    m.add_function(wrap_pyfunction!(line, &m)?)?;
    m.add_function(wrap_pyfunction!(brake, &m)?)?;
    m.add_function(wrap_pyfunction!(orbit_entry, &m)?)?;
    m.add_function(wrap_pyfunction!(orbit_exit, &m)?)?;
    m.add_function(wrap_pyfunction!(plane_through, &m)?)?;
    m.add_function(wrap_pyfunction!(rest, &m)?)?;
    m.add_function(wrap_pyfunction!(hyp_piece, &m)?)?;
    m.add_function(wrap_pyfunction!(lin_piece, &m)?)?;
    m.add_function(wrap_pyfunction!(null_piece, &m)?)?;
    m.add_function(wrap_pyfunction!(circ_piece, &m)?)?;
    m.setattr(
        "__all__",
        [
            "Brake",
            "End",
            "Line",
            "OrbitElements",
            "OrbitEntry",
            "Piece",
            "Plane",
            "Start",
            "brake",
            "circ_piece",
            "hyp_piece",
            "lin_piece",
            "line",
            "null_piece",
            "orbit_entry",
            "orbit_exit",
            "plane_through",
            "rest",
        ],
    )?;
    Ok(m)
}
