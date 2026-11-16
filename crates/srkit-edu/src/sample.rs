//! Turning a worldline into a polyline: which stretch of it is in view, and enough points along it.

use srkit_core::{kinematics::Event, roots::bisect, segment::Segment, worldline::Worldline};

use crate::{
    geometry::{Pt, Rect},
    view::Chart,
    Error,
};

/// The most points one worldline may take.
pub const MAX_POINTS: usize = 65_536;

/// The largest distance, in pixels, between a drawn worldline and the true curve.
pub const CURVE_TOLERANCE: f64 = 0.25;

/// A chord is accepted when the curve is within this many pixels of it at the quarter points and the middle.
/// Between those probes a smooth curve strays a little further, so the test is stricter than the promise.
const PROBE_TOLERANCE: f64 = 0.1;

/// Halvings allowed below the first cut of a segment.
const MAX_DEPTH: u32 = 40;

/// What a worldline is drawn onto: the chart, the window and the scale of each axis in pixels per unit.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Window {
    pub chart: Chart,
    pub rect: Rect,
    pub sx: f64,
    pub st: f64,
}

/// The stretch of coordinate time a worldline's segments cover: from the start of the first to the end of the
/// last. Outside it the worldline is only an analytic continuation, which is not drawn. A bound that is NaN
/// counts as open.
pub(crate) fn natural_range(worldline: &Worldline<'_>) -> (f64, f64) {
    let segments = worldline.segments();
    let lo = segments.first().map_or(f64::NAN, |a| a.seg.t0());
    let hi = segments.last().map_or(f64::NAN, |a| a.seg.t1());
    (
        if lo.is_nan() { f64::NEG_INFINITY } else { lo },
        if hi.is_nan() { f64::INFINITY } else { hi },
    )
}

fn event_at(worldline: &Worldline<'_>, t: f64) -> Event {
    Event::new(t, worldline.at(t).x)
}

/// The lab times `a..=b` between which the worldline is inside the window's time range, within `from..=to` and
/// the worldline's own range. `None` if no part is.
pub(crate) fn visible_range(
    worldline: &Worldline<'_>,
    window: &Window,
    from: f64,
    to: f64,
) -> Result<Option<(f64, f64)>, Error> {
    let (lo, hi) = natural_range(worldline);
    let (lo, hi) = (lo.max(from), hi.min(to));
    if !(lo <= hi) {
        return Ok(None);
    }
    let (t0, t1) = (window.rect.t0, window.rect.t1);
    let chart = &window.chart;
    if chart.rapidity == 0.0 {
        let (a, b) = (lo.max(t0), hi.min(t1));
        return Ok((a <= b).then_some((a, b)));
    }
    // In a boosted frame the time coordinate along the worldline is still an increasing function of lab time,
    // and between two of its events it changes by the lab time between them times a Doppler factor between
    // e^-|phi| and e^|phi|. That bounds the search.
    let time = |t: f64| chart.coords(event_at(worldline, t)).t;
    let start = chart.pivot.t.max(lo).min(hi);
    let at_start = time(start);
    let stretch = libm::exp(chart.rapidity.abs());
    let reach = (t0 - at_start).abs().max((t1 - at_start).abs()) * stretch * 1.000_001;
    if !at_start.is_finite() || !reach.is_finite() {
        return Err(Error::NonFinite);
    }
    let (a0, b0) = (lo.max(start - reach), hi.min(start + reach));
    let (at_a, at_b) = (time(a0), time(b0));
    if !at_a.is_finite() || !at_b.is_finite() {
        return Err(Error::NonFinite);
    }
    if at_b < t0 || at_a > t1 {
        return Ok(None);
    }
    let a = if at_a >= t0 {
        a0
    } else {
        bisect(|t| time(t) - t0, a0, b0, true)
    };
    let b = if at_b <= t1 {
        b0
    } else {
        bisect(|t| time(t) - t1, a0, b0, false)
    };
    Ok((a <= b).then_some((a, b)))
}

/// The worldline between lab times `a` and `b` as a polyline in chart coordinates.
///
/// Segments that are straight in spacetime (rest, inertial, light-like) take their two ends. Curved ones are
/// cut until the polyline is within [`CURVE_TOLERANCE`] pixels of the curve.
pub(crate) fn sample(
    worldline: &Worldline<'_>,
    window: &Window,
    a: f64,
    b: f64,
) -> Result<Vec<Pt>, Error> {
    let segments = worldline.segments();
    let mut sampler = Sampler {
        window,
        points: Vec::new(),
    };
    let (first, last) = (worldline.find(a), worldline.find(b));
    for i in first..=last {
        let Some(anchored) = segments.get(i) else {
            break;
        };
        let u = if i == first {
            a
        } else {
            anchored.seg.t0().max(a)
        };
        let v = match segments.get(i + 1) {
            Some(next) if i != last => next.seg.t0().min(b),
            _ => b,
        };
        if v >= u {
            sampler.segment(&anchored.seg, u, v)?;
        }
    }
    Ok(sampler.points)
}

struct Sampler<'a> {
    window: &'a Window,
    points: Vec<Pt>,
}

impl Sampler<'_> {
    fn eval(&self, segment: &Segment, t: f64) -> Result<Pt, Error> {
        let p = self.window.chart.coords(Event::new(t, segment.at(t).x));
        if p.is_finite() {
            Ok(p)
        } else {
            Err(Error::NonFinite)
        }
    }
    fn push(&mut self, p: Pt) -> Result<(), Error> {
        if self.points.last() == Some(&p) {
            return Ok(());
        }
        if self.points.len() >= MAX_POINTS {
            return Err(Error::TooDetailed);
        }
        self.points.push(p);
        Ok(())
    }
    /// How many equal pieces to cut a curved segment into before refining each. A circle gets at least eight
    /// to a turn, so that no bend can hide between the probes of one piece.
    fn pieces(segment: &Segment, u: f64, v: f64) -> Result<usize, Error> {
        let Segment::Circular(circle) = segment else {
            return Ok(4);
        };
        let eighths = libm::ceil((v - u) * circle.omega.abs() / core::f64::consts::FRAC_PI_4);
        if eighths.is_nan() {
            return Err(Error::NonFinite);
        }
        if eighths > MAX_POINTS as f64 {
            return Err(Error::TooDetailed);
        }
        Ok((eighths as usize).max(4))
    }
    fn segment(&mut self, segment: &Segment, u: f64, v: f64) -> Result<(), Error> {
        let (start, end) = (self.eval(segment, u)?, self.eval(segment, v)?);
        self.push(start)?;
        if matches!(
            segment,
            Segment::Rest(_) | Segment::Inertial(_) | Segment::Null(_)
        ) {
            return self.push(end);
        }
        let pieces = Self::pieces(segment, u, v)?;
        let (mut ta, mut pa) = (u, start);
        for j in 1..=pieces {
            let (tb, pb) = if j == pieces {
                (v, end)
            } else {
                let t = u + (v - u) * (j as f64 / pieces as f64);
                (t, self.eval(segment, t)?)
            };
            self.refine(segment, ta, pa, tb, pb, 0)?;
            (ta, pa) = (tb, pb);
        }
        Ok(())
    }
    /// Add the points after `pa` up to and including `pb`.
    fn refine(
        &mut self,
        segment: &Segment,
        ta: f64,
        pa: Pt,
        tb: f64,
        pb: Pt,
        depth: u32,
    ) -> Result<(), Error> {
        let tm = 0.5 * (ta + tb);
        if !(tm > ta && tm < tb) {
            // no representable time between the ends
            return self.push(pb);
        }
        let pm = self.eval(segment, tm)?;
        let probes = [
            self.eval(segment, 0.5 * (ta + tm))?,
            pm,
            self.eval(segment, 0.5 * (tm + tb))?,
        ];
        if probes.iter().all(|&q| self.near_chord(pa, pb, q)) {
            return self.push(pb);
        }
        if depth >= MAX_DEPTH {
            return Err(Error::TooDetailed);
        }
        self.refine(segment, ta, pa, tm, pm, depth + 1)?;
        self.refine(segment, tm, pm, tb, pb, depth + 1)
    }
    /// Whether `q` is within the probe tolerance of the chord from `a` to `b`, in pixels.
    fn near_chord(&self, a: Pt, b: Pt, q: Pt) -> bool {
        let (sx, st) = (self.window.sx, self.window.st);
        let (bx, by) = ((b.x - a.x) * sx, (b.t - a.t) * st);
        let (qx, qy) = ((q.x - a.x) * sx, (q.t - a.t) * st);
        let length2 = bx * bx + by * by;
        let along = if length2 > 0.0 {
            ((qx * bx + qy * by) / length2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let distance = libm::hypot(qx - along * bx, qy - along * by);
        // Far from the window the coordinates themselves are only known to a few units in the last place;
        // do not chase that noise.
        let size = (a.x.abs().max(b.x.abs()).max(q.x.abs())) * sx
            + (a.t.abs().max(b.t.abs()).max(q.t.abs())) * st;
        distance <= PROBE_TOLERANCE.max(16.0 * f64::EPSILON * size)
    }
}
