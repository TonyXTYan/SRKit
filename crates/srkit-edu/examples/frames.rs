//! One scene drawn in two frames: a station at rest and a ship passing it at 0.6 c.
//!
//! ```sh
//! cargo run -p srkit-edu --example frames > lab.svg
//! cargo run -p srkit-edu --example frames -- ship > ship.svg
//! cargo run -p srkit-edu --example frames -- ship dark > ship-dark.svg
//! ```
//!
//! In the lab frame the ship's axes are tilted towards the light cone. Drawn in the ship's frame they are
//! square and it is the station that moves; light is at 45° in both.

use std::process::ExitCode;

use srkit_core::{kinematics::Event, segment::Segment, worldline::WorldlineBuf, Vec3};
use srkit_edu::{
    Colour, Cone, Dash, Diagram, Error, Frame, Placement, Style, Theme, TickNumbers, View,
};

fn scene(from_ship: bool, theme: Theme) -> Result<String, Error> {
    let origin = Event::new(0.0, Vec3::ZERO);
    let velocity = Vec3::new(0.6, 0.0, 0.0);
    let station = WorldlineBuf::new([Segment::rest(-20.0, 20.0, Vec3::ZERO)], 0.0, 0.0)
        .map_err(|_| Error::NonFinite)?;
    let ship = WorldlineBuf::new(
        [Segment::inertial(-20.0, 20.0, velocity * -20.0, velocity)],
        0.0,
        0.0,
    )
    .map_err(|_| Error::NonFinite)?;
    let (station, ship) = (station.view(), ship.view());
    let ship_frame = Frame::moving_with(&ship, 0.0, Vec3::X)?;

    let mut view = View::equal_scales((-6.0, 6.0), (-2.0, 8.0), 48.0);
    if from_ship {
        view = view.frame(ship_frame);
    }
    let mut diagram = Diagram::new(view)?;
    diagram
        .theme(theme)
        .title(if from_ship {
            "The ship's frame: the station moves at 0.6 c"
        } else {
            "The station's frame: the ship moves at 0.6 c"
        })
        .axis_titles(
            if from_ship { "x′ (ly)" } else { "x (ly)" },
            if from_ship { "t′ (yr)" } else { "t (yr)" },
        );

    // What the ship sees of the station when the ship's clock reads 5 years.
    let now = ship.t_of_tau(5.0);
    let observer = Event::new(now, ship.at(now).x);
    diagram
        .frame_grid(
            ship_frame,
            2.0,
            "Ship's coordinate grid, every 2",
            Style::guide(),
        )?
        .frame_axes(
            ship_frame,
            "Ship's axes",
            Style::new(Colour::Series(1), 1.25, Dash::DashDot),
        )?
        .light_cone(
            origin,
            Cone::Both,
            "Light cone of the meeting",
            Style::light(),
        )?
        .seen(
            observer,
            &station,
            "Light by which the ship sees the station",
            Style::light().with_width(2.0),
        )?
        .worldline(&station, "Station", Style::series(0))?
        .worldline(&ship, "Ship", Style::series(1))?
        .proper_time_ticks(
            &station,
            1.0,
            if from_ship {
                TickNumbers::Left
            } else {
                TickNumbers::None
            },
            "Station clock, every year",
            Style::series(0),
        )?
        .proper_time_ticks(
            &ship,
            1.0,
            if from_ship {
                TickNumbers::None
            } else {
                TickNumbers::Right
            },
            "Ship clock, every year",
            Style::series(1),
        )?
        .event_placed(origin, "Meeting", Placement::Left, Style::horizon())?
        .event_placed(
            observer,
            "Ship looks back",
            Placement::Right,
            Style::horizon(),
        )?;
    Ok(diagram.render())
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let from_ship = arguments.iter().any(|a| a == "ship");
    let theme = if arguments.iter().any(|a| a == "dark") {
        Theme::Dark
    } else {
        Theme::Light
    };
    match scene(from_ship, theme) {
        Ok(svg) => {
            print!("{svg}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("frames: {error}");
            ExitCode::FAILURE
        }
    }
}
