//! Light-cone intersections on bounded, local-coordinate segments.
//!
//! A timelike source has a strictly decreasing retarded residual and a
//! timelike receiver has a strictly increasing future residual. Thus each
//! bounded query has at most one root; no simulation timestep is involved.
//! Null motion permits intervals of roots. Those use the latest retarded
//! event and the earliest future event within the segment.
//!
//! These are floating-point constructions. They preserve a numerical bracket
//! and return its causal side, but are **not** certified interval arithmetic or
//! exact integer arrival predicates. Rebase large epochs before querying.

use crate::kinematics::{Event, C};
use crate::segment::{Segment, State};
use crate::Error;

#[derive(Clone, Copy)]
enum Query {
    Retarded(Event),
    Future(Event),
}

impl Query {
    fn residual(self, state: &State) -> Result<f64, Error> {
        let value = match self {
            Self::Retarded(observer) => {
                (observer.t - state.event.t) - (observer.position - state.event.position).norm() / C
            }
            Self::Future(emission) => {
                (state.event.t - emission.t) - (state.event.position - emission.position).norm() / C
            }
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(Error::Overflow)
        }
    }
}

fn validate_event(event: Event) -> Result<(), Error> {
    if event.t.is_finite() && event.position.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidInput)
    }
}

/// Source state whose emitted light reaches `observer`.
///
/// The returned event lies within `segment` and no later than `observer`.
/// `NoIntersection` means the segment does not contain a retarded root; it does
/// not mean the source is invisible on every other segment of its worldline.
/// On a null interval coincident with the observer's past cone, chooses its
/// latest event. See the module-level floating-point limitations.
pub fn retarded_state(segment: &Segment, observer: Event) -> Result<State, Error> {
    validate_event(observer)?;
    let start = segment.start().t;
    let end = segment.end().t.min(observer.t);
    if end < start {
        return Err(Error::NoIntersection);
    }
    intersect(segment, Query::Retarded(observer), start, end)
}

/// Receiver state where light from `emission` first reaches `segment`.
///
/// The returned event lies within `segment` and no earlier than `emission`.
/// `NoIntersection` means no root lies inside this finite segment. On a null
/// interval coincident with the emission's future cone, chooses its earliest
/// event. See the module-level floating-point limitations.
pub fn future_intersection(segment: &Segment, emission: Event) -> Result<State, Error> {
    validate_event(emission)?;
    let start = segment.start().t.max(emission.t);
    let end = segment.end().t;
    if start > end {
        return Err(Error::NoIntersection);
    }
    intersect(segment, Query::Future(emission), start, end)
}

fn intersect(segment: &Segment, query: Query, mut lo: f64, mut hi: f64) -> Result<State, Error> {
    let left = query.residual(&segment.state(lo)?)?;
    let right = query.residual(&segment.state(hi)?)?;
    match query {
        Query::Retarded(_) => {
            if left < 0.0 || right > 0.0 {
                return Err(Error::NoIntersection);
            }
            if right == 0.0 {
                return segment.state(hi);
            }
        }
        Query::Future(_) => {
            if left > 0.0 || right < 0.0 {
                return Err(Error::NoIntersection);
            }
            if left == 0.0 {
                return segment.state(lo);
            }
        }
    }

    // 2200 halvings cover the finite f64 exponent range, including subnormals.
    // Usually adjacent representable times are reached after about 53 steps.
    for _ in 0..2200 {
        let difference = hi - lo;
        let mid = if difference.is_finite() {
            lo + difference * 0.5
        } else {
            lo * 0.5 + hi * 0.5
        };
        if mid <= lo || mid >= hi {
            return segment.state(match query {
                Query::Retarded(_) => lo,
                Query::Future(_) => hi,
            });
        }
        let residual = query.residual(&segment.state(mid)?)?;
        match query {
            Query::Retarded(_) if residual >= 0.0 => lo = mid,
            Query::Future(_) if residual < 0.0 => lo = mid,
            _ => hi = mid,
        }
    }
    Err(Error::PrecisionLoss)
}
