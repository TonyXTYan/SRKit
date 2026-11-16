//! Itineraries: sequences of segments that join continuously, with no instant change of velocity.
//!
//! A plan starts from a [`Start`] (a time, a proper time and a position) and is built from pieces, each a
//! closed-form [`Segment`] with its proper duration. [`line`] is the workhorse: a straight run that boosts,
//! coasts if the speed is capped, and brakes, between any two rapidities along the line. [`orbit_entry`] joins a
//! circular orbit along a tangent so that the velocity is continuous at the join.
//!
//! Nothing here allocates: a plan holds at most [`MAX_PIECES`] pieces.

use crate::{
    orbit::CircularOrbit,
    segment::{d_cosh, Circular, Hyperbolic, Inertial, Null, Rest, Segment, TauState},
    vec3::Vec3,
};

/// The most pieces a single plan function returns (an orbit entry from near the centre: boost, coast, brake,
/// boost along the tangent).
pub const MAX_PIECES: usize = 4;

/// What a piece is for. Applications label pieces from this.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    Rest,
    Boost,
    Coast,
    Brake,
    Orbit,
    Jump,
}

impl Phase {
    /// The phase as a lower-case word: `rest`, `boost`, `coast`, `brake`, `orbit`, `jump`.
    pub fn name(self) -> &'static str {
        match self {
            Self::Rest => "rest",
            Self::Boost => "boost",
            Self::Coast => "coast",
            Self::Brake => "brake",
            Self::Orbit => "orbit",
            Self::Jump => "jump",
        }
    }
}

/// Where and when a plan (or its next piece) begins.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Start {
    pub t: f64,
    pub tau: f64,
    pub x: Vec3,
}

/// The state at the end of a piece or a plan.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct End {
    pub t: f64,
    pub tau: f64,
    pub x: Vec3,
    pub v: Vec3,
    pub phi: f64,
    pub acc: f64,
}

impl End {
    /// Carry on from here.
    pub fn start(&self) -> Start {
        Start {
            t: self.t,
            tau: self.tau,
            x: self.x,
        }
    }
    fn at_rest(s: Start) -> Self {
        Self {
            t: s.t,
            tau: s.tau,
            x: s.x,
            v: Vec3::ZERO,
            phi: 0.0,
            acc: 0.0,
        }
    }
}

/// One piece of a plan: a segment, the proper time at its start, how much proper time it takes, and its purpose.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Piece {
    pub seg: Segment,
    pub tau0: f64,
    pub dtau: f64,
    pub phase: Phase,
}

impl Piece {
    /// The state at the end of the piece.
    pub fn end(&self) -> End {
        let e: TauState = self.seg.at_tau(self.dtau);
        End {
            t: e.t,
            tau: self.tau0 + self.dtau,
            x: e.x,
            v: e.v,
            phi: e.phi,
            acc: e.acc,
        }
    }
}

/// Up to [`MAX_PIECES`] pieces, in order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pieces {
    items: [Option<Piece>; MAX_PIECES],
    len: usize,
}

impl Default for Pieces {
    fn default() -> Self {
        Self {
            items: [None; MAX_PIECES],
            len: 0,
        }
    }
}

impl Pieces {
    fn push(&mut self, p: Piece) {
        self.items[self.len] = Some(p);
        self.len += 1;
    }
    fn extend(&mut self, other: &Pieces) {
        for p in other.iter() {
            self.push(*p);
        }
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    pub fn iter(&self) -> impl Iterator<Item = &Piece> {
        self.items[..self.len].iter().flatten()
    }
    pub fn get(&self, i: usize) -> Option<&Piece> {
        self.items.get(i).and_then(|p| p.as_ref())
    }
}

/// Stay at rest for proper time `dtau`.
pub fn rest(start: Start, dtau: f64) -> Piece {
    Piece {
        seg: Segment::Rest(Rest::new(start.t, start.t + dtau, start.x)),
        tau0: start.tau,
        dtau,
        phase: Phase::Rest,
    }
}

/// Accelerate at proper acceleration `alpha` (signed) along the unit vector `n`, from rapidity `phi0` to `phi1`
/// along `n`.
pub fn hyp_piece(start: Start, n: Vec3, alpha: f64, phi0: f64, phi1: f64, phase: Phase) -> Piece {
    let dtau = (phi1 - phi0) / alpha;
    let mut g = Hyperbolic::new(start.t, f64::NAN, start.x, n, alpha, phi0);
    let seg = Segment::Hyperbolic(g);
    g.t1 = seg.t_of_tau(dtau);
    Piece {
        seg: Segment::Hyperbolic(g),
        tau0: start.tau,
        dtau,
        phase,
    }
}

/// Coast at velocity `v` for proper time `dtau`.
pub fn lin_piece(start: Start, v: Vec3, dtau: f64) -> Piece {
    let mut g = Inertial::new(start.t, f64::NAN, start.x, v);
    g.t1 = start.t + dtau / g.ig;
    Piece {
        seg: Segment::Inertial(g),
        tau0: start.tau,
        dtau,
        phase: Phase::Coast,
    }
}

/// A light-like jump to `x1` (a teleport that moves at `c`): it takes the light-travel time in coordinate time
/// and no proper time.
pub fn null_piece(start: Start, x1: Vec3) -> Piece {
    Piece {
        seg: Segment::Null(Null::between(start.t, start.x, x1)),
        tau0: start.tau,
        dtau: 0.0,
        phase: Phase::Jump,
    }
}

/// The elements of a circle to fly: centre, radius, angular velocity, plane and the angle at the join.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitElements {
    pub centre: Vec3,
    pub r: f64,
    pub omega: f64,
    pub e1: Vec3,
    pub e2: Vec3,
    pub a0: f64,
    pub orbit: CircularOrbit,
}

/// Fly a circle for proper time `dtau` (which may be infinite), starting at its join angle.
pub fn circ_piece(start: Start, orb: &OrbitElements, dtau: f64) -> Piece {
    let mut g = Circular::new(
        start.t,
        f64::NAN,
        orb.centre,
        orb.r,
        orb.omega,
        orb.e1,
        orb.e2,
        orb.a0,
    );
    g.t1 = if dtau.is_finite() {
        start.t + dtau / g.ig
    } else {
        f64::INFINITY
    };
    Piece {
        seg: Segment::Circular(g),
        tau0: start.tau,
        dtau,
        phase: Phase::Orbit,
    }
}

/// A straight run and its summary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line {
    pub pieces: Pieces,
    pub end: End,
    /// Unit vector along the run.
    pub n: Vec3,
    /// Total change of rapidity the drive supplies (boost plus brake): what the rocket equation charges for.
    pub dphi: f64,
    /// Length of the run.
    pub length: f64,
    /// Peak rapidity.
    pub phi_peak: f64,
    /// Length coasted at the cap (zero if the cap is not reached).
    pub coast_length: f64,
}

/// `acosh(1 + x)` for `x >= 0`, accurate for small `x`: `2 asinh(sqrt(x / 2))`.
pub fn acosh1p(x: f64) -> f64 {
    2.0 * libm::asinh(libm::sqrt(0.5 * x.max(0.0)))
}

/// `cosh(phi) - 1` without cancellation: `2 sinh²(phi / 2)`.
pub fn cosh_m1(phi: f64) -> f64 {
    let s = libm::sinh(0.5 * phi);
    2.0 * s * s
}

/// A straight run from `start.x`, where the rapidity along the line is `phi0 >= 0`, to `x1`, where it is
/// `phi1 >= 0`: boost at proper acceleration `a`, coast if the rapidity reaches `phi_cap`, brake at `a_brake`.
///
/// Returns `None` if the run is too short to change speed as asked. `n_default` is the direction to use when the
/// run has no length but a change of speed is asked for.
///
/// With no coast the peak rapidity `phi_m` satisfies
/// `(cosh phi_m - cosh phi0) / a + (cosh phi_m - cosh phi1) / a_brake = L`. It is solved for `cosh phi_m - 1`, so
/// a hop of a few light-seconds keeps as many digits as a crossing of the galaxy.
#[allow(clippy::too_many_arguments)]
pub fn line(
    start: Start,
    x1: Vec3,
    a: f64,
    phi_cap: f64,
    phi0: f64,
    phi1: f64,
    a_brake: f64,
    n_default: Option<Vec3>,
) -> Option<Line> {
    let d = x1 - start.x;
    let length = d.norm();
    if length < 1e-12 && phi0 == 0.0 && phi1 == 0.0 {
        return Some(Line {
            pieces: Pieces::default(),
            end: End::at_rest(start),
            n: n_default.unwrap_or(Vec3::X),
            dphi: 0.0,
            length: 0.0,
            phi_peak: 0.0,
            coast_length: 0.0,
        });
    }
    let n = if length > 0.0 {
        d / length
    } else {
        n_default.unwrap_or(Vec3::X)
    };
    let b2 = a_brake;
    // cosh(phi_m) - 1
    let chm1 = (length + cosh_m1(phi0) / a + cosh_m1(phi1) / b2) / (1.0 / a + 1.0 / b2);
    let mut phm = acosh1p(chm1);
    let mut coast_length = 0.0;
    if phm > phi_cap {
        phm = phi_cap;
        coast_length =
            length - d_cosh(phi0, phi_cap - phi0) / a - d_cosh(phi1, phi_cap - phi1) / b2;
    }
    if phm < phi0.max(phi1) - 1e-12 {
        return None;
    }
    let mut pieces = Pieces::default();
    let mut s = start;
    let mut end = End {
        t: start.t,
        tau: start.tau,
        x: start.x,
        v: n * libm::tanh(phi0),
        phi: phi0,
        acc: 0.0,
    };
    if phm > phi0 + 1e-15 {
        let p = hyp_piece(s, n, a, phi0, phm, Phase::Boost);
        pieces.push(p);
        end = p.end();
        s = end.start();
    }
    if coast_length > 1e-12 {
        let b = libm::tanh(phm);
        let p = lin_piece(s, n * b, coast_length / b * libm::sqrt(1.0 - b * b));
        pieces.push(p);
        end = p.end();
        s = end.start();
    }
    if phm > phi1 + 1e-15 {
        let p = hyp_piece(s, n, -b2, phm, phi1, Phase::Brake);
        pieces.push(p);
        end = p.end();
    }
    // Pin the end on the target: rounding in the closed forms is far below a metre at these scales.
    end.x = x1;
    Some(Line {
        pieces,
        end,
        n,
        dphi: (phm - phi0) + (phm - phi1),
        length,
        phi_peak: phm,
        coast_length,
    })
}

/// Brake from rapidity `phi0` along the unit vector `n` to rest, at proper acceleration `a`.
/// Returns the piece (none if already at rest), the end state and the rapidity spent.
pub fn brake(start: Start, n: Vec3, phi0: f64, a: f64) -> (Option<Piece>, End, f64) {
    if !(phi0 > 0.0) {
        return (None, End::at_rest(start), 0.0);
    }
    let p = hyp_piece(start, n, -a, phi0, 0.0, Phase::Brake);
    (Some(p), p.end(), phi0)
}

/// An orthonormal pair `(e1, e2)` and their normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    pub e1: Vec3,
    pub e2: Vec3,
    pub normal: Vec3,
}

fn norm_or_self(v: Vec3) -> Vec3 {
    let n = v.norm();
    if n > 0.0 {
        v / n
    } else {
        v
    }
}

/// The plane through `centre` that holds the point `x` (or the direction `e1_default` when `x` is at the centre)
/// and is tilted least from the `xy` plane. `e1` points from the centre to `x`.
pub fn plane_through(centre: Vec3, x: Vec3, e1_default: Option<Vec3>) -> Plane {
    let d = x - centre;
    let dist = d.norm();
    if dist < 1e-12 {
        let e1 = norm_or_self(e1_default.unwrap_or(Vec3::X));
        let mut normal = norm_or_self(Vec3::Z - e1 * e1.z);
        if e1.cross(Vec3::Z).norm() < 1e-9 {
            normal = e1.perpendicular();
        }
        return Plane {
            e1,
            e2: normal.cross(e1),
            normal,
        };
    }
    let e1 = d / dist;
    let raw = Vec3::Z - e1 * e1.z;
    let normal = if raw.norm() < 1e-9 {
        e1.perpendicular()
    } else {
        norm_or_self(raw)
    };
    Plane {
        e1,
        e2: normal.cross(e1),
        normal,
    }
}

/// Options for [`orbit_entry`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EntryOptions {
    /// The direction of `e1` when the start is at the centre.
    pub e1: Option<Vec3>,
    /// The join angle when entering from near the centre (default 0).
    pub a0: Option<f64>,
    /// Proper acceleration of the boost along the tangent, when it differs from the hop's.
    pub a_tangent: Option<f64>,
}

/// The pieces that bring a ship onto a circular orbit, and the circle's elements at the join.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitEntry {
    pub pieces: Pieces,
    pub end: End,
    pub orb: OrbitElements,
    pub dphi: f64,
    /// True when the entry started near the centre (a hop to a point behind the join, then a boost along the
    /// tangent); false for the tangent run from outside.
    pub from_inside: bool,
}

/// Enter a circular orbit of radius `r` at speed `beta` about `centre`, from rest at `start.x`, at proper
/// acceleration `a`. The orbit runs counter-clockwise in `(e1, e2)`.
///
/// From outside, the run is the straight line tangent to the circle at its join (boost, coast, brake to the
/// orbit's speed), so the velocity is continuous. From near the centre: a rest-to-rest run to a point behind the
/// join, then a boost along the tangent. The circle itself is added by the caller with [`circ_piece`], so that it
/// can be open-ended. `None` if a run is too short for its change of speed.
pub fn orbit_entry(
    start: Start,
    centre: Vec3,
    r: f64,
    beta: f64,
    a: f64,
    phi_cap: f64,
    opts: EntryOptions,
) -> Option<OrbitEntry> {
    let ob = CircularOrbit::from_radius_speed(r, beta);
    let pl = plane_through(centre, start.x, opts.e1);
    let dist = start.x.dist(centre);
    let a_t = match opts.a_tangent {
        Some(v) if v != 0.0 && !v.is_nan() => v,
        _ => a,
    };
    let db = cosh_m1(ob.phi) / a_t;
    let mut pieces = Pieces::default();
    let mut dphi = 0.0;
    if dist > r * 1.000001 && libm::sqrt(dist * dist - r * r) >= db * 1.05 {
        // Outside: the tangent from the start, which sits at (dist, 0) in the plane; the join is at angle theta
        // with cos theta = r / dist.
        let th = libm::acos(r / dist);
        let q = centre + pl.e1 * (r * libm::cos(th)) + pl.e2 * (r * libm::sin(th));
        let l = line(start, q, a, phi_cap.max(ob.phi), 0.0, ob.phi, a, None)?;
        pieces.extend(&l.pieces);
        dphi += l.dphi;
        Some(OrbitEntry {
            pieces,
            end: l.end,
            orb: OrbitElements {
                centre,
                r,
                omega: ob.omega,
                e1: pl.e1,
                e2: pl.e2,
                a0: th,
                orbit: ob,
            },
            dphi,
            from_inside: false,
        })
    } else {
        // Near the centre: the join q at angle a0, and a point s behind it on the tangent.
        let a0 = opts.a0.unwrap_or(0.0);
        let q = centre + pl.e1 * (r * libm::cos(a0)) + pl.e2 * (r * libm::sin(a0));
        let tangent = pl.e2 * libm::cos(a0) - pl.e1 * libm::sin(a0);
        let s_point = q - tangent * db;
        let l1 = line(start, s_point, a, phi_cap, 0.0, 0.0, a, None)?;
        pieces.extend(&l1.pieces);
        dphi += l1.dphi;
        let p = hyp_piece(l1.end.start(), tangent, a_t, 0.0, ob.phi, Phase::Boost);
        pieces.push(p);
        dphi += ob.phi;
        let mut end = p.end();
        end.x = q;
        Some(OrbitEntry {
            pieces,
            end,
            orb: OrbitElements {
                centre,
                r,
                omega: ob.omega,
                e1: pl.e1,
                e2: pl.e2,
                a0,
                orbit: ob,
            },
            dphi,
            from_inside: true,
        })
    }
}

/// Leave a circle moving at velocity `v`: brake along the tangent to rest at proper acceleration `a`.
pub fn orbit_exit(start: Start, v: Vec3, a: f64) -> (Option<Piece>, End, f64) {
    let b = v.norm();
    if !(b > 0.0) {
        // already at rest: nothing to brake
        return brake(start, Vec3::X, 0.0, a);
    }
    brake(start, v / b, libm::atanh(b), a)
}
