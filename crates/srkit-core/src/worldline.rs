//! Worldlines: ordered segments with a proper-time anchor on each.
//!
//! The first segment is taken to run back to the infinite past and the last on to the infinite future, so a
//! worldline answers every query; give it a closing [`crate::segment::Rest`] if the body should stay put. Give
//! the first segment a finite `t0`: it is the reference its proper time is counted from, and queries before it
//! use the segment's analytic continuation.
//!
//! [`Worldline`] borrows its segments and allocates nothing. With the `alloc` feature, [`WorldlineBuf`] owns
//! them.

use crate::{
    kinematics::Event,
    segment::{Segment, State, TauState},
    Error,
};

/// A segment and the body's proper time at its start.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Anchored {
    pub seg: Segment,
    pub tau0: f64,
}

impl Anchored {
    pub const fn new(seg: Segment, tau0: f64) -> Self {
        Self { seg, tau0 }
    }
}

/// Set every segment's `tau0` so that proper time is continuous along the chain and equals `anchor_tau` at
/// coordinate time `anchor_t`.
pub fn anchor(segs: &mut [Anchored], anchor_t: f64, anchor_tau: f64) {
    if segs.is_empty() {
        return;
    }
    let k = find_index(segs, anchor_t);
    segs[k].tau0 = anchor_tau - segs[k].seg.dtau(anchor_t);
    for i in k + 1..segs.len() {
        segs[i].tau0 = segs[i - 1].tau0 + segs[i - 1].seg.dtau(segs[i].seg.t0());
    }
    for i in (0..k).rev() {
        segs[i].tau0 = segs[i + 1].tau0 - segs[i].seg.dtau(segs[i + 1].seg.t0());
    }
}

/// Index of the last segment starting at or before `t` (0 if `t` is before them all).
fn find_index(segs: &[Anchored], t: f64) -> usize {
    last_le(segs.len(), t, |i| segs[i].seg.t0())
}

/// Index of the last of `n` ascending keys that is `<= x`, or 0 if none is.
fn last_le(n: usize, x: f64, key: impl Fn(usize) -> f64) -> usize {
    if n == 0 || x < key(0) {
        return 0;
    }
    let (mut lo, mut hi) = (0, n - 1);
    while lo < hi {
        let m = (lo + hi).div_ceil(2);
        if key(m) <= x {
            lo = m;
        } else {
            hi = m - 1;
        }
    }
    lo
}

/// A worldline over borrowed, anchored segments.
#[derive(Clone, Copy, Debug)]
pub struct Worldline<'a> {
    segs: &'a [Anchored],
    has_null: bool,
}

impl<'a> Worldline<'a> {
    /// Wrap anchored segments (see [`anchor`]). They must be in order of `t0`; an empty slice is an error.
    pub fn new(segs: &'a [Anchored]) -> Result<Self, Error> {
        if segs.is_empty() {
            return Err(Error::InvalidInput);
        }
        if segs.windows(2).any(|w| !(w[0].seg.t0() <= w[1].seg.t0())) {
            return Err(Error::InvalidInput);
        }
        Ok(Self {
            segs,
            has_null: segs.iter().any(|a| a.seg.is_null()),
        })
    }
    pub fn segments(&self) -> &'a [Anchored] {
        self.segs
    }
    pub fn len(&self) -> usize {
        self.segs.len()
    }
    pub fn is_empty(&self) -> bool {
        self.segs.is_empty()
    }
    /// Whether any segment is light-like (a jump).
    pub fn has_null(&self) -> bool {
        self.has_null
    }
    /// Index of the segment in force at coordinate time `t`.
    pub fn find(&self, t: f64) -> usize {
        find_index(self.segs, t)
    }
    /// Index of the segment in force at proper time `tau`.
    pub fn find_tau(&self, tau: f64) -> usize {
        last_le(self.segs.len(), tau, |i| self.segs[i].tau0)
    }
    /// Position and velocity at coordinate time `t`.
    pub fn at(&self, t: f64) -> State {
        self.segs[self.find(t)].seg.at(t)
    }
    /// Proper time at coordinate time `t`.
    pub fn tau(&self, t: f64) -> f64 {
        let a = &self.segs[self.find(t)];
        a.tau0 + a.seg.dtau(t)
    }
    /// The state at proper time `tau`.
    pub fn at_tau(&self, tau: f64) -> TauState {
        let a = &self.segs[self.find_tau(tau)];
        a.seg.at_tau(tau - a.tau0)
    }
    /// Coordinate time at proper time `tau`.
    pub fn t_of_tau(&self, tau: f64) -> f64 {
        let a = &self.segs[self.find_tau(tau)];
        a.seg.t_of_tau(tau - a.tau0)
    }
    /// The largest speed on the worldline.
    pub fn vmax(&self) -> f64 {
        self.segs.iter().fold(0.0, |m, a| m.max(a.seg.vmax()))
    }

    fn end_of(&self, i: usize) -> f64 {
        if i + 1 < self.segs.len() {
            self.segs[i + 1].seg.t0()
        } else {
            f64::INFINITY
        }
    }
    fn start_of(&self, i: usize) -> f64 {
        if i == 0 {
            f64::NEG_INFINITY
        } else {
            self.segs[i].seg.t0()
        }
    }

    /// The coordinate time at which light from `event` reaches this worldline, or `+inf` if it never does (a
    /// body that outruns the light, such as one accelerating away for ever).
    ///
    /// `f(t) = t - event.t - |x(t) - event.x|` increases along the worldline, so a binary search over segment
    /// starts finds the segment and [`Segment::recv`] finds the root within it. Where light can run alongside a
    /// light-like segment `f` is flat, and the search takes the first root.
    pub fn recv(&self, event: Event) -> f64 {
        let n = self.segs.len();
        let f = |i: usize| {
            let g = &self.segs[i].seg;
            g.t0() - event.t - g.at(g.t0()).x.dist(event.x)
        };
        let (mut a, mut b) = (self.find(event.t), n - 1);
        while a < b {
            let m = (a + b).div_ceil(2);
            let inside = if self.has_null {
                f(m) < 0.0
            } else {
                f(m) <= 0.0
            };
            if inside {
                a = m;
            } else {
                b = m - 1;
            }
        }
        for i in a..n {
            if let Some(t) = self.segs[i]
                .seg
                .recv(event, self.start_of(i), self.end_of(i))
            {
                return t;
            }
        }
        f64::INFINITY
    }

    /// The coordinate time of the event on this worldline whose light reaches `event`: what an observer there
    /// sees of this body. `-inf` if no such moment exists.
    pub fn send(&self, event: Event) -> f64 {
        let f = |i: usize| {
            let g = &self.segs[i].seg;
            g.t0() + event.x.dist(g.at(g.t0()).x) - event.t
        };
        let (mut a, mut b) = (0, self.find(event.t));
        while a < b {
            let m = (a + b).div_ceil(2);
            if f(m) <= 0.0 {
                a = m;
            } else {
                b = m - 1;
            }
        }
        for i in (0..=a).rev() {
            if let Some(t) = self.segs[i]
                .seg
                .send(event, self.start_of(i), self.end_of(i))
            {
                return t;
            }
        }
        f64::NEG_INFINITY
    }

    /// The state this body is seen in from `event` (its position and velocity when the light left), with the
    /// emission time, or `None` if nothing of it can be seen.
    pub fn seen_from(&self, event: Event) -> Option<(f64, State)> {
        let t = self.send(event);
        t.is_finite().then(|| (t, self.at(t)))
    }

    /// Check that consecutive segments join: each starts where and when the one before ends, at the same
    /// velocity. A light-like segment is exempt from the velocity check on both sides (a jump starts and ends at
    /// rest by convention).
    pub fn check_continuity(&self, position_tol: f64, velocity_tol: f64) -> Result<(), Error> {
        for w in self.segs.windows(2) {
            let (a, b) = (&w[0].seg, &w[1].seg);
            let t = b.t0();
            let (end, start) = (a.at(t), b.at(t));
            let end_x = if let Segment::Null(g) = a {
                g.p1
            } else {
                end.x
            };
            if !(end_x.dist(start.x) <= position_tol) {
                return Err(Error::Discontinuous);
            }
            if !a.is_null() && !b.is_null() && !(end.v.dist(start.v) <= velocity_tol) {
                return Err(Error::Discontinuous);
            }
            let tau_end = w[0].tau0 + a.dtau(t);
            if !((tau_end - w[1].tau0).abs() <= position_tol.max(1e-9 * tau_end.abs())) {
                return Err(Error::Discontinuous);
            }
        }
        Ok(())
    }
}

#[cfg(feature = "alloc")]
mod owned {
    use super::{anchor, Anchored, Worldline};
    use crate::{segment::Segment, Error};
    use alloc::vec::Vec;

    /// A worldline that owns its segments.
    #[derive(Clone, Debug, Default, PartialEq)]
    pub struct WorldlineBuf {
        segs: Vec<Anchored>,
    }

    impl WorldlineBuf {
        /// Build from segments in order of `t0`, anchoring proper time to `anchor_tau` at `anchor_t`.
        pub fn new(
            segments: impl IntoIterator<Item = Segment>,
            anchor_t: f64,
            anchor_tau: f64,
        ) -> Result<Self, Error> {
            let mut segs: Vec<Anchored> = segments
                .into_iter()
                .map(|seg| Anchored::new(seg, 0.0))
                .collect();
            if segs.is_empty() {
                return Err(Error::InvalidInput);
            }
            anchor(&mut segs, anchor_t, anchor_tau);
            Worldline::new(&segs)?;
            Ok(Self { segs })
        }
        /// Build from segments that already carry their proper-time anchors.
        pub fn from_anchored(segs: Vec<Anchored>) -> Result<Self, Error> {
            Worldline::new(&segs)?;
            Ok(Self { segs })
        }
        /// Append a segment that starts where the last one ends in coordinate time; its proper-time anchor is
        /// computed from the chain.
        pub fn push(&mut self, seg: Segment) -> Result<(), Error> {
            let last = self.segs.last().ok_or(Error::InvalidInput)?;
            if !(seg.t0() >= last.seg.t0()) {
                return Err(Error::InvalidInput);
            }
            let tau0 = last.tau0 + last.seg.dtau(seg.t0());
            self.segs.push(Anchored::new(seg, tau0));
            Ok(())
        }
        /// Borrow as a [`Worldline`] for queries.
        pub fn view(&self) -> Worldline<'_> {
            Worldline::new(&self.segs).expect("validated when built")
        }
        pub fn segments(&self) -> &[Anchored] {
            &self.segs
        }
    }
}

#[cfg(feature = "alloc")]
pub use owned::WorldlineBuf;
