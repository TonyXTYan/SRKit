//! Two ready-made diagrams, built with the same public calls a caller would use.
//!
//! Both work in years and light-years, so an acceleration is in light-years per year squared: one standard
//! gravity is 1.0323.
//!
//! Special relativity has no scale of its own: a trip of `d` at acceleration `a` is the trip of 1 at `a d`,
//! magnified `d` times. Each scene is worked out in such a unit trip and scaled up at the end, because the
//! core's light-cone solver stops at a few parts in 10¹⁶ of one unit and would be coarse on a scene whose
//! numbers are all tiny.

use srkit_core::{
    kinematics::Event,
    plan::{self, Start},
    segment::Segment,
    worldline::WorldlineBuf,
    Vec3,
};

use crate::{
    axis::Step,
    diagram::{Diagram, Placement, TickNumbers},
    style::{Colour, Dash, Style},
    svg::fixed,
    view::View,
    Error,
};

/// A scene parameter must be a positive, finite number.
fn positive(value: f64) -> Result<(), Error> {
    if !value.is_finite() {
        Err(Error::NonFinite)
    } else if value > 0.0 {
        Ok(())
    } else {
        Err(Error::OutOfRange)
    }
}

/// A number for a title. Numbers of ordinary size are written with `decimals` decimals; very large and very
/// small ones as a power of ten, `1.00e-9`, so that a title cannot run to hundreds of digits.
fn readable(value: f64, decimals: usize) -> String {
    if value == 0.0 || (1e-3..1e6).contains(&value.abs()) {
        fixed(value, decimals)
    } else {
        format!("{value:.2e}")
    }
}

/// As [`readable`] with up to four decimals and no trailing zeros: `4.2465`, `1.5`, `2`.
fn trimmed(value: f64) -> String {
    let text = readable(value, 4);
    if text.contains('e') || !text.contains('.') {
        return text;
    }
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

fn worldline(segments: impl IntoIterator<Item = Segment>) -> Result<WorldlineBuf, Error> {
    WorldlineBuf::new(segments, 0.0, 0.0).map_err(|_| Error::NonFinite)
}

/// The same segment with every time and length multiplied by `k`. Rapidities and directions are unchanged
/// and proper accelerations are divided by `k`.
fn magnified(segment: &Segment, k: f64) -> Segment {
    match segment {
        Segment::Rest(g) => Segment::rest(g.t0 * k, g.t1 * k, g.p * k),
        Segment::Inertial(g) => Segment::inertial(g.t0 * k, g.t1 * k, g.p0 * k, g.v),
        Segment::Hyperbolic(g) => {
            Segment::hyperbolic(g.t0 * k, g.t1 * k, g.p0 * k, g.n, g.alpha / k, g.phi0)
        }
        // A straight run with no speed cap is made of the three kinds above.
        other => *other,
    }
}

fn magnified_event(event: Event, k: f64) -> Event {
    Event::new(event.t * k, event.x * k)
}

/// The twin paradox as an SVG diagram.
///
/// One twin stays at home. The other flies `distance` light-years out and back at a constant proper
/// acceleration of `acceleration` light-years per year squared, turning the ship round at each half-way point
/// so that it stops at the far end and again at home. Both clocks are ticked, and each twin sends the other a
/// flash of light at every tick of their own clock, so the diagram shows what each sees of the other's ageing.
///
/// `twin_paradox(4.2465, 1.0323)` is Proxima Centauri and back at one gravity. Light runs at 45° unless the
/// trip is so slow that the plot would be a sliver, in which case the position axis is stretched.
///
/// Both numbers must be positive and finite. [`Error::OutOfRange`] also covers a trip too slow to tell from
/// Newton's (`distance × acceleration` below about 10⁻³⁰).
pub fn twin_paradox(distance: f64, acceleration: f64) -> Result<String, Error> {
    Ok(twin_paradox_diagram(distance, acceleration)?.render())
}

/// The [`twin_paradox`] diagram before it is rendered, for a caller who wants another theme or more on it.
pub fn twin_paradox_diagram(distance: f64, acceleration: f64) -> Result<Diagram, Error> {
    positive(distance)?;
    positive(acceleration)?;
    // The unit trip: one unit out and one back, at proper acceleration `a d`.
    let unit_acceleration = acceleration * distance;
    positive(unit_acceleration)?;
    let run = |from: Start, to: Vec3| {
        plan::line(
            from,
            to,
            unit_acceleration,
            f64::INFINITY,
            0.0,
            0.0,
            unit_acceleration,
            None,
        )
        .ok_or(Error::OutOfRange)
    };
    let launch = Start {
        t: 0.0,
        tau: 0.0,
        x: Vec3::ZERO,
    };
    let out = run(launch, Vec3::X)?;
    let back = run(out.end.start(), Vec3::ZERO)?;
    if out.pieces.is_empty() || back.pieces.is_empty() {
        return Err(Error::OutOfRange);
    }
    let unit_segments: Vec<Segment> = out
        .pieces
        .iter()
        .chain(back.pieces.iter())
        .map(|piece| piece.seg)
        .collect();
    let unit_ship = worldline(unit_segments.iter().copied())?;
    let unit_ship = unit_ship.view();

    // The real trip.
    let (t_turn, t_end, tau_end) = (
        out.end.t * distance,
        back.end.t * distance,
        back.end.tau * distance,
    );
    if !(t_end.is_finite() && tau_end.is_finite()) {
        return Err(Error::NonFinite);
    }
    if !(t_end > 0.0 && tau_end > 0.0) {
        return Err(Error::OutOfRange);
    }
    let far = Vec3::new(distance, 0.0, 0.0);
    let ship = worldline(unit_segments.iter().map(|s| magnified(s, distance)))?;
    let home = worldline([Segment::rest(0.0, t_end, Vec3::ZERO)])?;
    let (ship, home) = (ship.view(), home.view());

    // Equal scales, unless that would leave the plot too narrow to read. There is room on the right for the
    // ship's clock readings and the label of the turnaround.
    let x = (-0.2 * distance, 1.5 * distance);
    let t = (-0.04 * t_end, 1.04 * t_end);
    let height = 640.0;
    let width = ((x.1 - x.0) * (height / (t.1 - t.0))).max(280.0);
    let mut diagram = Diagram::new(View::new(x, t, width, height))?;
    diagram
        .title(&format!(
            "Twin paradox: {} ly out and back at {} ly/yr²",
            trimmed(distance),
            trimmed(acceleration)
        ))
        .subtitle(&format!(
            "At the reunion the home clock reads {} yr and the ship clock {} yr",
            readable(t_end, 2),
            readable(tau_end, 2)
        ))
        .axis_titles("x (ly)", "t (yr)");

    // Between six and twelve ticks on each clock.
    let home_step = Step::at_least(t_end / 12.0);
    let ship_step = Step::at_least(tau_end / 12.0);
    let (home_style, ship_style) = (Style::series(0), Style::series(1));
    let home_light = Style::light().with_colour(home_style.colour);
    let ship_light = Style::light()
        .with_colour(ship_style.colour)
        .with_dash(Dash::Dotted);

    diagram
        .worldline(&home, "Twin at home", home_style)?
        .worldline(&ship, "Travelling twin", ship_style)?
        // The home clock keeps coordinate time, so the time axis already numbers its ticks.
        .proper_time_ticks(
            &home,
            home_step.value(),
            TickNumbers::None,
            &format!("Home clock, every {} yr", home_step.text()),
            home_style,
        )?
        .proper_time_ticks(
            &ship,
            ship_step.value(),
            TickNumbers::Right,
            &format!("Ship clock, every {} yr", ship_step.text()),
            ship_style,
        )?;
    for k in 1..=12 {
        let sent = f64::from(k) * home_step.value();
        if !(sent < t_end * (1.0 - 1e-9)) {
            break;
        }
        // When the ship receives it, found on the unit trip.
        let received = unit_ship.recv(Event::new(sent / distance, Vec3::ZERO));
        if received <= back.end.t {
            let arrival = Event::new(received, unit_ship.at(received).x);
            diagram.signal(
                Event::new(sent, Vec3::ZERO),
                magnified_event(arrival, distance),
                "Light from home, one flash per tick",
                home_light,
            )?;
        }
    }
    for k in 1..=12 {
        let tau = f64::from(k) * ship_step.value();
        if !(tau < tau_end * (1.0 - 1e-9)) {
            break;
        }
        // Home is at the origin, so the flash arrives after the light-travel time from where the ship is.
        let sent = ship.t_of_tau(tau);
        let from = Event::new(sent, ship.at(sent).x);
        let to = Event::new(sent + from.x.norm(), Vec3::ZERO);
        diagram.signal(
            from,
            to,
            "Light from the ship, one flash per tick",
            ship_light,
        )?;
    }
    diagram
        .event_placed(
            Event::new(t_turn, far),
            "Turnaround",
            Placement::Right,
            Style::horizon(),
        )?
        .event_placed(
            Event::new(t_end, Vec3::ZERO),
            "Reunion",
            Placement::Above,
            Style::horizon(),
        )?;
    Ok(diagram)
}

/// A ship that accelerates for ever, the horizon that follows it, and light from home that does and does not
/// catch it, as an SVG diagram.
///
/// The ship leaves home from rest at `t = 0` and holds a proper acceleration of `acceleration` light-years
/// per year squared. Its Rindler horizon is the light ray that trails it by `1 / acceleration` light-years:
/// a flash sent from home before `t = 1 / acceleration` years reaches the ship, and one sent at or after that
/// moment never does while the engine runs.
///
/// The acceleration must be positive and finite, and between 10⁻¹⁰⁰ and 10¹⁰⁰.
pub fn rindler_chase(acceleration: f64) -> Result<String, Error> {
    Ok(rindler_chase_diagram(acceleration)?.render())
}

/// The [`rindler_chase`] diagram before it is rendered.
pub fn rindler_chase_diagram(acceleration: f64) -> Result<Diagram, Error> {
    positive(acceleration)?;
    // The horizon distance, and the last moment a flash from home can still reach the ship.
    let h = 1.0 / acceleration;
    if !(1e-100..=1e100).contains(&h) {
        return Err(Error::OutOfRange);
    }
    let x = (-1.5 * h, 3.7 * h);
    let t = (-0.5 * h, 4.3 * h);
    let mut diagram = Diagram::new(View::equal_scales(x, t, 560.0 / (t.1 - t.0)))?;
    diagram
        .title(&format!(
            "Rindler chase: a ship holding {} ly/yr²",
            trimmed(acceleration)
        ))
        .subtitle(&format!(
            "Light sent from home at or after t = {} yr never reaches the ship",
            readable(h, 2)
        ))
        .axis_titles("x (ly)", "t (yr)");

    let drive = Segment::hyperbolic(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, acceleration, 0.0);
    let ship = worldline([Segment::rest(t.0, 0.0, Vec3::ZERO), drive])?;
    let home = worldline([Segment::rest(t.0, f64::INFINITY, Vec3::ZERO)])?;
    let (ship, home) = (ship.view(), home.view());
    // The unit ship: proper acceleration 1, horizon one unit behind it.
    let unit_ship = worldline([Segment::hyperbolic(
        0.0,
        f64::INFINITY,
        Vec3::ZERO,
        Vec3::X,
        1.0,
        0.0,
    )])?;
    let unit_ship = unit_ship.view();

    let step = Step::at_least(libm::asinh(t.1 * acceleration) * h / 5.0);
    diagram
        .worldline(&home, "Home", Style::series(0))?
        .worldline(&ship, "Ship", Style::series(1))?
        .proper_time_ticks(
            &ship,
            step.value(),
            TickNumbers::Right,
            &format!("Ship clock, every {} yr", step.text()),
            Style::series(1),
        )?
        .rindler_horizon(&drive, "Rindler horizon of the ship", Style::horizon())?;

    let caught = Style::light();
    let lost = Style::light()
        .with_colour(Colour::Muted)
        .with_dash(Dash::Dotted);
    for fraction in [0.3, 0.6, 0.8, 1.2, 1.6] {
        let from = Event::new(fraction * h, Vec3::ZERO);
        let received = unit_ship.recv(Event::new(fraction, Vec3::ZERO));
        if received.is_finite() {
            let to = magnified_event(Event::new(received, unit_ship.at(received).x), h);
            diagram
                .signal(from, to, "Light from home that reaches the ship", caught)?
                .event(to, "", caught)?;
        } else {
            // It trails the ship for ever: draw it to a point far outside the view.
            let far = 16.0 * h;
            let to = Event::new(from.t + far, Vec3::new(far, 0.0, 0.0));
            diagram.signal(from, to, "Light from home that never arrives", lost)?;
        }
    }
    diagram.event_placed(
        Event::new(h, Vec3::ZERO),
        "Last chance",
        Placement::Left,
        Style::horizon(),
    )?;
    Ok(diagram)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_are_trimmed() {
        assert_eq!(trimmed(4.2465), "4.2465");
        assert_eq!(trimmed(1.5), "1.5");
        assert_eq!(trimmed(2.0), "2");
        assert_eq!(trimmed(100.0), "100");
        assert_eq!(trimmed(0.001), "0.001");
        assert_eq!(trimmed(0.00001), "1.00e-5");
        assert_eq!(trimmed(2.5e9), "2.50e9");
        assert_eq!(readable(11.7454, 2), "11.75");
        assert_eq!(readable(1e90, 2), "1.00e90");
        assert_eq!(readable(0.0, 2), "0.00");
    }

    #[test]
    fn a_magnified_segment_is_the_same_motion_at_another_scale() {
        let unit = Segment::hyperbolic(0.5, 3.0, Vec3::new(0.2, 0.0, 0.0), Vec3::X, -0.7, 1.1);
        let big = magnified(&unit, 40.0);
        for t in [0.5, 1.0, 2.25, 3.0] {
            let (a, b) = (unit.at(t), big.at(40.0 * t));
            assert!((a.x * 40.0).dist(b.x) < 1e-12 && a.v.dist(b.v) < 1e-14);
            assert!((unit.dtau(t) * 40.0 - big.dtau(40.0 * t)).abs() < 1e-12);
        }
        let rest = magnified(&Segment::rest(1.0, 2.0, Vec3::new(3.0, 0.0, 0.0)), 0.5);
        assert_eq!(rest, Segment::rest(0.5, 1.0, Vec3::new(1.5, 0.0, 0.0)));
    }
}
