//! The window a diagram looks through: ranges, plot size, the spatial axis and the frame it is drawn in.

use srkit_core::{
    kinematics::{boost_rapidity, Event},
    worldline::Worldline,
    Vec3,
};

use crate::{geometry::Pt, Error};

/// The largest plot side, in pixels.
pub const MAX_PIXELS: f64 = 16384.0;

/// The largest view bound, in either direction.
const MAX_BOUND: f64 = 1e150;

/// An inertial frame moving along the diagram's axis: an origin event and a rapidity.
///
/// The frame's coordinates agree with the lab's at the origin event and are boosted about it, so a frame of
/// rapidity zero is the lab whatever its origin. With the origin at the lab's own origin this is the textbook
/// Lorentz transformation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    /// The event the boost leaves in place.
    pub origin: Event,
    /// Rapidity of the frame along the axis direction: `atanh` of its velocity as a fraction of `c`. Negative
    /// for a frame moving against the axis.
    pub rapidity: f64,
}

impl Frame {
    /// A frame boosted by `rapidity` along the axis about `origin`.
    pub const fn new(origin: Event, rapidity: f64) -> Self {
        Self { origin, rapidity }
    }
    /// A frame moving at `velocity` (a fraction of `c`, signed along the axis) about `origin`.
    /// [`Error::OutOfRange`] at or above the speed of light.
    pub fn from_velocity(origin: Event, velocity: f64) -> Result<Self, Error> {
        if velocity.is_nan() {
            return Err(Error::NonFinite);
        }
        if !(velocity.abs() < 1.0) {
            return Err(Error::OutOfRange);
        }
        Ok(Self::new(origin, libm::atanh(velocity)))
    }
    /// The rest frame of the body on `worldline` at coordinate time `t`: its origin is the body's event then,
    /// and its rapidity is that of the body's velocity along `direction`, the diagram's axis.
    pub fn moving_with(worldline: &Worldline<'_>, t: f64, direction: Vec3) -> Result<Self, Error> {
        let direction = unit(direction)?;
        let state = worldline.at(t);
        if !t.is_finite() || !state.x.is_finite() || !state.v.is_finite() {
            return Err(Error::NonFinite);
        }
        Self::from_velocity(Event::new(t, state.x), state.v.dot(direction))
    }
}

/// The unit vector along `direction`.
fn unit(direction: Vec3) -> Result<Vec3, Error> {
    if !direction.is_finite() {
        return Err(Error::NonFinite);
    }
    direction.normalized().map_err(|_| Error::ZeroDirection)
}

/// What a diagram shows and how large it is drawn.
///
/// The ranges are in the coordinates of the frame the diagram is drawn in: position along the axis across, time
/// up. `width` and `height` are the size of the plot area in pixels; the axes, title and legend are added
/// around it. Use [`View::equal_scales`] for a diagram in which light runs at 45°.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct View {
    /// Left and right ends of the position range.
    pub x: (f64, f64),
    /// Bottom and top ends of the time range.
    pub t: (f64, f64),
    /// Width of the plot area in pixels.
    pub width: f64,
    /// Height of the plot area in pixels.
    pub height: f64,
    /// The point in space where the position coordinate is zero.
    pub origin: Vec3,
    /// The direction in space of increasing position. It is scaled to unit length.
    pub direction: Vec3,
    /// The frame the diagram is drawn in. `None` is the lab frame, in which the worldlines are stored.
    pub frame: Option<Frame>,
}

impl View {
    /// A view of positions `x.0..x.1` and times `t.0..t.1` on a plot of `width` by `height` pixels, along the
    /// lab's `x` axis from its origin, in the lab frame.
    pub const fn new(x: (f64, f64), t: (f64, f64), width: f64, height: f64) -> Self {
        Self {
            x,
            t,
            width,
            height,
            origin: Vec3::ZERO,
            direction: Vec3::X,
            frame: None,
        }
    }
    /// A view with the same scale on both axes, `pixels_per_unit` pixels to one unit of time and to the
    /// distance light covers in it, so that light runs at 45°.
    pub fn equal_scales(x: (f64, f64), t: (f64, f64), pixels_per_unit: f64) -> Self {
        Self::new(
            x,
            t,
            (x.1 - x.0) * pixels_per_unit,
            (t.1 - t.0) * pixels_per_unit,
        )
    }
    /// The same view along another spatial axis: positions are measured from `origin` along `direction`.
    /// Anything off the axis is projected onto it.
    pub const fn axis(mut self, origin: Vec3, direction: Vec3) -> Self {
        self.origin = origin;
        self.direction = direction;
        self
    }
    /// The same view drawn in `frame`. The ranges then mean that frame's coordinates.
    pub const fn frame(mut self, frame: Frame) -> Self {
        self.frame = Some(frame);
        self
    }
    /// Check the view and scale its direction to unit length.
    pub(crate) fn validated(mut self) -> Result<Self, Error> {
        let numbers = [
            self.x.0,
            self.x.1,
            self.t.0,
            self.t.1,
            self.width,
            self.height,
        ];
        if numbers.iter().any(|v| !v.is_finite()) || !self.origin.is_finite() {
            return Err(Error::NonFinite);
        }
        self.direction = unit(self.direction)?;
        if let Some(frame) = self.frame {
            check_frame(frame)?;
        }
        for (lo, hi) in [self.x, self.t] {
            let size = lo.abs().max(hi.abs());
            if size > MAX_BOUND {
                return Err(Error::TooLarge);
            }
            // A range must be wide enough for its ticks to be told apart in binary64.
            if !(hi - lo > 1e-9 * size) || !(hi - lo > 1.0 / MAX_BOUND) {
                return Err(Error::EmptyView);
            }
        }
        for side in [self.width, self.height] {
            if side < 1.0 {
                return Err(Error::EmptyView);
            }
            if side > MAX_PIXELS {
                return Err(Error::TooLarge);
            }
        }
        Ok(self)
    }
}

pub(crate) fn check_frame(frame: Frame) -> Result<(), Error> {
    if frame.origin.is_finite() && frame.rapidity.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFinite)
    }
}

/// The map from lab events to the coordinates of a frame along the axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Chart {
    /// Where the position coordinate is zero.
    pub origin: Vec3,
    /// Unit vector along the axis.
    pub dir: Vec3,
    /// The frame's fixed event.
    pub pivot: Event,
    /// The fixed event's position coordinate.
    pub pivot_x: f64,
    /// Rapidity of the frame along `dir`.
    pub rapidity: f64,
}

impl Chart {
    /// The chart of `frame` (the lab if `None`) along the axis from `origin` in the unit direction `dir`.
    pub fn new(origin: Vec3, dir: Vec3, frame: Option<Frame>) -> Self {
        let frame = frame.unwrap_or(Frame::new(Event::new(0.0, origin), 0.0));
        Self {
            origin,
            dir,
            pivot: frame.origin,
            pivot_x: (frame.origin.x - origin).dot(dir),
            rapidity: frame.rapidity,
        }
    }
    /// The frame's coordinates of a lab event.
    pub fn coords(&self, event: Event) -> Pt {
        if self.rapidity == 0.0 {
            return Pt {
                x: (event.x - self.origin).dot(self.dir),
                t: event.t,
            };
        }
        let (t, x) = boost_rapidity(
            self.dir * self.rapidity,
            event.t - self.pivot.t,
            event.x - self.pivot.x,
        );
        Pt {
            x: self.pivot_x + x.dot(self.dir),
            t: self.pivot.t + t,
        }
    }
    /// The lab event on the axis with the given coordinates in the frame.
    pub fn event(&self, p: Pt) -> Event {
        if self.rapidity == 0.0 {
            return Event::new(p.t, self.origin + self.dir * p.x);
        }
        let (t, x) = boost_rapidity(
            self.dir * -self.rapidity,
            p.t - self.pivot.t,
            self.dir * (p.x - self.pivot_x),
        );
        Event::new(self.pivot.t + t, self.origin + self.dir * self.pivot_x + x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_chart_and_its_inverse_agree() {
        let dir = Vec3::new(0.6, 0.8, 0.0);
        let origin = Vec3::new(1.0, -2.0, 0.5);
        let frame = Frame::new(Event::new(3.0, origin + dir * 2.0), 0.7);
        let chart = Chart::new(origin, dir, Some(frame));
        // the fixed event keeps its lab coordinates
        let fixed = chart.coords(frame.origin);
        assert!((fixed.x - 2.0).abs() < 1e-15 && fixed.t == 3.0);
        for (x, t) in [(0.0, 0.0), (5.0, -3.0), (-7.5, 11.0)] {
            let back = chart.coords(chart.event(Pt { x, t }));
            assert!((back.x - x).abs() < 1e-12 && (back.t - t).abs() < 1e-12);
        }
        // a body at rest in the frame keeps its position coordinate
        let v = dir * 0.7_f64.tanh();
        let (a, b) = (
            chart.coords(Event::new(0.0, origin)),
            chart.coords(Event::new(4.0, origin + v * 4.0)),
        );
        assert!((a.x - b.x).abs() < 1e-12 && b.t > a.t);
    }

    #[test]
    fn rapidity_zero_is_the_lab_whatever_the_origin() {
        let lab = Chart::new(Vec3::ZERO, Vec3::X, None);
        let same = Chart::new(
            Vec3::ZERO,
            Vec3::X,
            Some(Frame::new(Event::new(9.0, Vec3::new(4.0, 1.0, 0.0)), 0.0)),
        );
        let event = Event::new(1.25, Vec3::new(-3.5, 2.0, 0.0));
        assert_eq!(lab.coords(event), same.coords(event));
        assert_eq!(lab.coords(event), Pt { x: -3.5, t: 1.25 });
    }

    #[test]
    fn views_are_checked() {
        let good = View::new((0.0, 1.0), (0.0, 1.0), 100.0, 100.0);
        assert!(good.validated().is_ok());
        let scaled = View::equal_scales((-1.0, 3.0), (0.0, 8.0), 50.0);
        assert_eq!((scaled.width, scaled.height), (200.0, 400.0));
        let direction = good
            .axis(Vec3::ZERO, Vec3::new(0.0, 3.0, 4.0))
            .validated()
            .unwrap()
            .direction;
        assert!((direction.norm() - 1.0).abs() < 1e-15);
        let cases = [
            (
                View::new((1.0, 1.0), (0.0, 1.0), 100.0, 100.0),
                Error::EmptyView,
            ),
            (
                View::new((2.0, 1.0), (0.0, 1.0), 100.0, 100.0),
                Error::EmptyView,
            ),
            (
                View::new((0.0, 1.0), (1e12, 1e12 + 1e-6), 100.0, 100.0),
                Error::EmptyView,
            ),
            (
                View::new((0.0, 1.0), (0.0, 1.0), 0.5, 100.0),
                Error::EmptyView,
            ),
            (
                View::new((0.0, 1.0), (0.0, 1.0), -100.0, 100.0),
                Error::EmptyView,
            ),
            (
                View::new((0.0, f64::NAN), (0.0, 1.0), 100.0, 100.0),
                Error::NonFinite,
            ),
            (
                View::new((0.0, 1.0), (0.0, f64::INFINITY), 100.0, 100.0),
                Error::NonFinite,
            ),
            (
                View::new((0.0, 1.0), (0.0, 1.0), f64::NAN, 100.0),
                Error::NonFinite,
            ),
            (
                View::new((0.0, 1e200), (0.0, 1.0), 100.0, 100.0),
                Error::TooLarge,
            ),
            (
                View::new((0.0, 1.0), (0.0, 1.0), 1e6, 100.0),
                Error::TooLarge,
            ),
            (good.axis(Vec3::ZERO, Vec3::ZERO), Error::ZeroDirection),
            (
                good.axis(Vec3::new(f64::NAN, 0.0, 0.0), Vec3::X),
                Error::NonFinite,
            ),
            (
                good.axis(Vec3::ZERO, Vec3::new(f64::INFINITY, 0.0, 0.0)),
                Error::NonFinite,
            ),
            (
                good.frame(Frame::new(Event::new(0.0, Vec3::ZERO), f64::NAN)),
                Error::NonFinite,
            ),
        ];
        for (view, error) in cases {
            assert_eq!(view.validated(), Err(error), "{view:?}");
        }
    }

    #[test]
    fn frames_from_velocities() {
        let origin = Event::new(0.0, Vec3::ZERO);
        let frame = Frame::from_velocity(origin, -0.6).unwrap();
        assert!((frame.rapidity + 0.6_f64.atanh()).abs() < 1e-15);
        assert_eq!(Frame::from_velocity(origin, 1.0), Err(Error::OutOfRange));
        assert_eq!(
            Frame::from_velocity(origin, f64::NAN),
            Err(Error::NonFinite)
        );
    }
}
