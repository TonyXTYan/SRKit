//! Bounded analytic segments in a local inertial coordinate frame (SI units).
//!
//! Hyperbolic motion uses the plane spanned by the initial four-velocity and
//! a unit thrust vector fixed in the initial comoving frame. This supports
//! arbitrary 3D thrust angles, not just collinear motion.
use crate::kinematics::{lorentz_factor, Event, Vec3, C};
use crate::Error;

#[derive(Clone, Copy, Debug)]
pub struct State {
    pub event: Event,
    pub velocity: Vec3,
    /// Elapsed proper seconds since the segment's start.
    pub proper_time: f64,
}

#[derive(Clone, Copy, Debug)]
enum Motion {
    Inertial {
        velocity: Vec3,
        gamma: f64,
    },
    Accelerated(Acceleration),
    Circular {
        center: Vec3,
        x: Vec3,
        y: Vec3,
        omega: f64,
        gamma: f64,
    },
    Null {
        direction: Vec3,
    },
}

#[derive(Clone, Copy, Debug)]
struct Acceleration {
    initial_velocity: Vec3,
    gamma: f64,
    speed: f64,
    axis: Vec3,
    cosine: f64,
    perpendicular: Vec3,
    // Positive null components of the temporal boost. Avoid subtracting
    // nearly equal gamma-sized quantities for relativistic counter-thrust.
    plus: f64,
    minus: f64,
    a: f64,
}

impl Acceleration {
    fn new(initial_velocity: Vec3, gamma: f64, direction: Vec3, a: f64) -> Self {
        let beta = initial_velocity / C;
        let speed = beta.norm();
        let (axis, cosine, perpendicular) = if speed == 0.0 {
            (direction, 1.0, Vec3::ZERO)
        } else {
            let axis = beta / speed;
            let cosine = axis.dot(direction);
            let perpendicular = direction - axis * cosine;
            // Renormalize in this basis so roundoff near collinearity cannot
            // create a direction cosine outside [-1,1].
            let length = libm::hypot(cosine, perpendicular.norm());
            (axis, cosine / length, perpendicular / length)
        };
        let (plus, minus) = if cosine >= 0.0 {
            (1.0 + speed * cosine, (1.0 - speed) + speed * (1.0 - cosine))
        } else {
            ((1.0 - speed) + speed * (1.0 + cosine), 1.0 - speed * cosine)
        };
        Self {
            initial_velocity,
            gamma,
            speed,
            axis,
            cosine,
            perpendicular,
            plus,
            minus,
            a,
        }
    }

    /// Temporal, parallel and perpendicular displacement divided by rapidity.
    /// The limit at zero keeps tiny burns accurate without forming C²/a.
    fn displacement_ratios(self, q: f64) -> (f64, f64, f64) {
        let half = 0.5 * q;
        let sinhc_half = if half == 0.0 {
            1.0
        } else {
            libm::sinh(half) / half
        };
        let cm1_over_q = half * sinhc_half * sinhc_half;
        let (temporal, parallel) = if q < 0.5 {
            let sinhc = if q == 0.0 { 1.0 } else { libm::sinh(q) / q };
            (
                sinhc + self.speed * self.cosine * cm1_over_q,
                self.speed * sinhc + self.cosine * cm1_over_q,
            )
        } else {
            let positive = libm::expm1(q) / q;
            let negative = -libm::expm1(-q) / q;
            (
                0.5 * (self.plus * positive + self.minus * negative),
                0.5 * ((self.speed + self.cosine) * positive
                    + (self.speed - self.cosine) * negative),
            )
        };
        (temporal, parallel, cm1_over_q)
    }

    fn elapsed(self, proper_time: f64) -> f64 {
        let q = self.a * proper_time / C;
        proper_time * (self.gamma * self.displacement_ratios(q).0)
    }

    fn rapidity(self, dt: f64) -> f64 {
        let s = (self.a * dt / C) / self.gamma;
        let bn = self.speed * self.cosine;
        let scale = libm::sqrt(self.plus * self.minus);
        if s < 0.5 {
            let root = libm::hypot(s + bn, scale);
            libm::log1p(2.0 * s / (1.0 - s + root))
        } else {
            let shift = 0.5 * (libm::log(self.plus) - libm::log(self.minus));
            libm::asinh((s + bn) / scale) - shift
        }
    }

    fn evaluate(self, dt: f64) -> (Vec3, Vec3, f64) {
        if dt == 0.0 {
            return (Vec3::ZERO, self.initial_velocity, 0.0);
        }
        let q = self.rapidity(dt);
        let (temporal, parallel, perpendicular) = self.displacement_ratios(q);
        // dt = gamma * tau * temporal. Unlike C/a*q this also remains finite
        // when acceleration is so small that C/a overflows.
        let tau = dt / (self.gamma * temporal);
        let displacement =
            (self.axis * (self.gamma * parallel) + self.perpendicular * perpendicular) * (C * tau);
        let (denominator, numerator, sinh) = if q < 0.5 {
            let sinh = libm::sinh(q);
            let cosh = libm::cosh(q);
            (
                cosh + self.speed * self.cosine * sinh,
                self.speed * cosh + self.cosine * sinh,
                sinh,
            )
        } else {
            let positive = libm::exp(q);
            let negative = libm::exp(-q);
            (
                0.5 * (self.plus * positive + self.minus * negative),
                0.5 * ((self.speed + self.cosine) * positive
                    + (self.speed - self.cosine) * negative),
                0.5 * (positive - negative),
            )
        };
        let velocity = (self.axis * (numerator / denominator)
            + self.perpendicular * (sinh / (self.gamma * denominator)))
            * C;
        (displacement, velocity, tau)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Segment {
    origin: Event,
    duration: f64,
    proper_duration: f64,
    motion: Motion,
}

fn valid_interval(origin: Event, duration: f64) -> Result<(), Error> {
    if !origin.t.is_finite()
        || !origin.position.is_finite()
        || !duration.is_finite()
        || duration <= 0.0
    {
        return Err(Error::InvalidInput);
    }
    if !(origin.t + duration).is_finite() {
        return Err(Error::Overflow);
    }
    if origin.t + duration <= origin.t {
        return Err(Error::PrecisionLoss);
    }
    Ok(())
}

impl Segment {
    pub fn inertial(origin: Event, velocity: Vec3, duration: f64) -> Result<Self, Error> {
        valid_interval(origin, duration)?;
        let gamma = lorentz_factor(velocity)?;
        let out = Self {
            origin,
            duration,
            proper_duration: duration / gamma,
            motion: Motion::Inertial { velocity, gamma },
        };
        out.check_end()
    }

    /// Constant positive proper acceleration in a fixed comoving direction.
    /// Duration is proper seconds; use `inertial` for zero acceleration.
    ///
    /// Positive null-coordinate components preserve counter-thrust accuracy at
    /// high initial gamma. Local floating-point time still limits the inverse
    /// near a high-gamma turnaround. A velocity rounding to c is rejected with
    /// `PrecisionLoss`, rather than clamped to an invented velocity.
    pub fn accelerated(
        origin: Event,
        initial_velocity: Vec3,
        rest_direction: Vec3,
        proper_acceleration: f64,
        proper_duration: f64,
    ) -> Result<Self, Error> {
        if !proper_acceleration.is_finite()
            || proper_acceleration <= 0.0
            || !proper_duration.is_finite()
            || proper_duration <= 0.0
        {
            return Err(Error::InvalidInput);
        }
        if proper_acceleration * proper_duration / C == 0.0 {
            // A zero rounded rapidity can still hide a representable
            // displacement; do not silently replace the burn by inertial motion.
            return Err(Error::PrecisionLoss);
        }
        let gamma = lorentz_factor(initial_velocity)?;
        let n = rest_direction.normalized()?;
        let acceleration = Acceleration::new(initial_velocity, gamma, n, proper_acceleration);
        let duration = acceleration.elapsed(proper_duration);
        valid_interval(origin, duration)?;
        let out = Self {
            origin,
            duration,
            proper_duration,
            motion: Motion::Accelerated(acceleration),
        };
        out.check_end()
    }

    /// A laboratory-frame circular orbit. `axis_x` and `axis_y` are orthogonal
    /// radius vectors in metres. The origin must equal center + axis_x.
    /// Phase is local floating point; use exact modular clocks for deep epochs.
    pub fn circular(
        origin: Event,
        center: Vec3,
        axis_x: Vec3,
        axis_y: Vec3,
        angular_velocity: f64,
        duration: f64,
    ) -> Result<Self, Error> {
        valid_interval(origin, duration)?;
        let r = axis_x.norm();
        if !center.is_finite()
            || !axis_x.is_finite()
            || !axis_y.is_finite()
            || !r.is_finite()
            || !axis_y.norm().is_finite()
            || !angular_velocity.is_finite()
            || r <= 0.0
            || (axis_y.norm() / r - 1.0).abs() > 1e-12
            || (axis_x / r).dot(axis_y / r).abs() > 1e-12
            || (origin.position - center - axis_x).norm() > 1e-12 * r
        {
            return Err(Error::InvalidInput);
        }
        let gamma = lorentz_factor(axis_y * angular_velocity)?;
        if !(angular_velocity * duration).is_finite() {
            return Err(Error::Overflow);
        }
        let out = Self {
            origin,
            duration,
            proper_duration: duration / gamma,
            motion: Motion::Circular {
                center,
                x: axis_x,
                y: axis_y,
                omega: angular_velocity,
                gamma,
            },
        };
        out.check_end()
    }

    pub fn null(origin: Event, direction: Vec3, duration: f64) -> Result<Self, Error> {
        valid_interval(origin, duration)?;
        let direction = direction.normalized()?;
        let out = Self {
            origin,
            duration,
            proper_duration: 0.0,
            motion: Motion::Null { direction },
        };
        out.check_end()
    }

    fn check_end(self) -> Result<Self, Error> {
        self.state(self.origin.t + self.duration)?;
        Ok(self)
    }

    pub fn start(&self) -> Event {
        self.origin
    }
    pub fn end(&self) -> Event {
        // Constructors validate this evaluation.
        self.state(self.origin.t + self.duration)
            .expect("validated segment")
            .event
    }
    pub fn duration(&self) -> f64 {
        self.duration
    }
    pub fn proper_duration(&self) -> f64 {
        self.proper_duration
    }
    pub fn is_null(&self) -> bool {
        matches!(self.motion, Motion::Null { .. })
    }

    /// Invert elapsed proper time to coordinate time. A null clock is not
    /// invertible; even tau=0 has a continuum of coordinate times.
    pub fn coord_time(&self, proper_time: f64) -> Result<f64, Error> {
        if self.is_null()
            || !proper_time.is_finite()
            || proper_time < 0.0
            || proper_time > self.proper_duration
        {
            return Err(Error::OutOfRange);
        }
        let dt = match self.motion {
            Motion::Accelerated(acceleration) => acceleration.elapsed(proper_time),
            _ => self.duration * (proper_time / self.proper_duration),
        };
        let t = self.origin.t + dt;
        if t.is_finite() {
            Ok(t)
        } else {
            Err(Error::Overflow)
        }
    }

    pub fn proper_time(&self, t: f64) -> Result<f64, Error> {
        Ok(self.state(t)?.proper_time)
    }

    pub fn state(&self, t: f64) -> Result<State, Error> {
        if !t.is_finite() || t < self.origin.t || t > self.origin.t + self.duration {
            return Err(Error::OutOfRange);
        }
        let dt = t - self.origin.t;
        let (position, velocity, proper_time) = match self.motion {
            Motion::Inertial { velocity, gamma } => {
                (self.origin.position + velocity * dt, velocity, dt / gamma)
            }
            Motion::Null { direction } => (
                self.origin.position + direction * (C * dt),
                direction * C,
                0.0,
            ),
            Motion::Circular {
                center,
                x,
                y,
                omega,
                gamma,
            } => {
                let phase = omega * dt;
                let (s, c) = (libm::sin(phase), libm::cos(phase));
                (center + x * c + y * s, (y * c - x * s) * omega, dt / gamma)
            }
            Motion::Accelerated(acceleration) => {
                let (displacement, velocity, tau) = acceleration.evaluate(dt);
                (self.origin.position + displacement, velocity, tau)
            }
        };
        if !position.is_finite() || !velocity.is_finite() || !proper_time.is_finite() {
            return Err(Error::Overflow);
        }
        if !self.is_null() && velocity.norm() >= C {
            return Err(Error::PrecisionLoss);
        }
        Ok(State {
            event: Event { t, position },
            velocity,
            proper_time,
        })
    }
}

/// A borrowed, allocation-free chain. Clock time resets per segment; this
/// wrapper sums proper durations. Joins must preserve position and velocity.
pub struct Worldline<'a> {
    segments: &'a [Segment],
}

impl<'a> Worldline<'a> {
    pub fn new(
        segments: &'a [Segment],
        position_tolerance: f64,
        velocity_tolerance: f64,
    ) -> Result<Self, Error> {
        if segments.is_empty()
            || !position_tolerance.is_finite()
            || position_tolerance < 0.0
            || !velocity_tolerance.is_finite()
            || velocity_tolerance < 0.0
        {
            return Err(Error::InvalidInput);
        }
        for pair in segments.windows(2) {
            let end = pair[0].state(pair[0].end().t)?;
            let start = pair[1].state(pair[1].start().t)?;
            if end.event.t != start.event.t
                || (end.event.position - start.event.position).norm() > position_tolerance
                || (end.velocity - start.velocity).norm() > velocity_tolerance
            {
                return Err(Error::Discontinuous);
            }
        }
        Ok(Self { segments })
    }

    pub fn state(&self, t: f64) -> Result<State, Error> {
        let mut elapsed = 0.0;
        for segment in self.segments {
            if t >= segment.start().t && t <= segment.end().t {
                let mut state = segment.state(t)?;
                state.proper_time += elapsed;
                return Ok(state);
            }
            elapsed += segment.proper_duration();
        }
        Err(Error::OutOfRange)
    }

    pub fn coord_time(&self, tau: f64) -> Result<f64, Error> {
        if !tau.is_finite() || tau < 0.0 {
            return Err(Error::OutOfRange);
        }
        let mut elapsed = 0.0;
        for segment in self.segments {
            if tau <= elapsed + segment.proper_duration() {
                return segment.coord_time(tau - elapsed);
            }
            elapsed += segment.proper_duration();
        }
        Err(Error::OutOfRange)
    }
}
