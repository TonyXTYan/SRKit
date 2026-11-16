//! Geometry: the coordinates written into the SVG are read back and checked against the physics.
mod common;

use common::{at_rest, distance_to_polylines, find, moving, parse, single, Element};
use srkit_core::{
    kinematics::Event,
    rindler::RindlerHorizon,
    segment::{Burn, Hyperbolic, Segment, Steered},
    worldline::{Worldline, WorldlineBuf},
    Vec3,
};
use srkit_edu::{Cone, Diagram, Error, Frame, Layout, Style, TickNumbers, View, CURVE_TOLERANCE};

/// Coordinates are written with two decimals, so a difference of two of them is good to a hundredth.
const ROUNDING: f64 = 0.011;

fn origin() -> Event {
    Event::new(0.0, Vec3::ZERO)
}

fn on_axis(t: f64, x: f64) -> Event {
    Event::new(t, Vec3::new(x, 0.0, 0.0))
}

fn square_view() -> View {
    View::equal_scales((-6.0, 6.0), (-4.0, 8.0), 40.0)
}

fn elements(diagram: &Diagram) -> Vec<Element> {
    parse(&diagram.render()).unwrap()
}

/// `(dx, dy)` of a line in pixels.
fn run(line: &Element) -> (f64, f64) {
    let ((x1, y1), (x2, y2)) = line.ends();
    (x2 - x1, y2 - y1)
}

fn close(a: (f64, f64), b: (f64, f64), tolerance: f64) -> bool {
    (a.0 - b.0).abs() <= tolerance && (a.1 - b.1).abs() <= tolerance
}

/// The pixel position of a lab event.
fn pixel(diagram: &Diagram, event: Event) -> (f64, f64) {
    let (x, t) = diagram.coordinates(event);
    diagram.layout().to_pixels(x, t)
}

/// The coordinates of a pixel position in the drawing frame.
fn coordinates(layout: &Layout, p: (f64, f64)) -> (f64, f64) {
    (
        layout.x.0 + (p.0 - layout.plot_left) / layout.plot_width * (layout.x.1 - layout.x.0),
        layout.t.1 - (p.1 - layout.plot_top) / layout.plot_height * (layout.t.1 - layout.t.0),
    )
}

#[test]
fn light_cones_are_at_45_degrees_with_equal_scales() {
    let frames = [
        None,
        Some(Frame::new(origin(), 0.9)),
        Some(Frame::new(on_axis(1.0, 2.0), -1.7)),
    ];
    for frame in frames {
        let mut view = square_view();
        view.frame = frame;
        let mut inside = 0;
        let apexes = [
            (0.0, 0.0),
            (2.5, -1.25),
            (-1.0, 2.0),
            (1.0, 0.5),
            (4.0, 3.0),
            (3.0, -2.0),
            (1.0, 2.0),
            (1.5, 2.2),
            (0.8, 1.7),
        ];
        for (t, x) in apexes {
            let apex = on_axis(t, x);
            let mut diagram = Diagram::new(view).unwrap();
            diagram
                .light_cone(apex, Cone::Both, "light", Style::light())
                .unwrap();
            // an apex outside the view shows only the rays that cross it
            let (ax, at) = diagram.coordinates(apex);
            if !(ax > -5.5 && ax < 5.5 && at > -3.5 && at < 7.5) {
                continue;
            }
            inside += 1;
            let all = elements(&diagram);
            let lines = find(&all, "cone", "line");
            assert_eq!(lines.len(), 4, "{frame:?} {apex:?}");
            let at = pixel(&diagram, apex);
            let mut directions = Vec::new();
            for line in lines {
                let (dx, dy) = run(line);
                assert!(dx.abs() > 15.0, "a ray of {dx} px");
                assert!(
                    (dx.abs() - dy.abs()).abs() <= 2.0 * ROUNDING,
                    "slope {} in {frame:?}",
                    dy / dx
                );
                // every ray starts at the apex
                assert!(close(line.ends().0, at, ROUNDING));
                directions.push((dx > 0.0, dy > 0.0));
            }
            directions.sort();
            directions.dedup();
            assert_eq!(directions.len(), 4, "one ray into each quadrant");
        }
        assert!(inside >= 3, "{frame:?}: only {inside} apexes in view");
    }
    // only the future, or only the past
    for (cone, up) in [(Cone::Future, true), (Cone::Past, false)] {
        let mut diagram = Diagram::new(square_view()).unwrap();
        diagram
            .light_cone(origin(), cone, "", Style::light())
            .unwrap();
        let all = elements(&diagram);
        let lines = find(&all, "cone", "line");
        assert_eq!(lines.len(), 2);
        // pixel y grows downwards
        assert!(lines.iter().all(|l| (run(l).1 < 0.0) == up));
    }
}

#[test]
fn with_unequal_scales_light_follows_the_scales() {
    // 30 px to a unit of position, 60 px to a unit of time
    let mut diagram = Diagram::new(View::new((-5.0, 5.0), (0.0, 5.0), 300.0, 300.0)).unwrap();
    diagram
        .light_cone(on_axis(1.0, 0.0), Cone::Future, "", Style::light())
        .unwrap();
    let all = elements(&diagram);
    for line in find(&all, "cone", "line") {
        let (dx, dy) = run(line);
        assert!((dy.abs() / dx.abs() - 2.0).abs() < 1e-3);
    }
}

#[test]
fn a_body_at_rest_is_vertical_and_a_moving_one_leans_by_its_speed() {
    let mut diagram = Diagram::new(square_view()).unwrap();
    let (rest, ship) = (at_rest(1.5, -100.0, 100.0), moving(0.6, -100.0, 100.0));
    diagram
        .worldline(&rest.view(), "rest", Style::series(0))
        .unwrap()
        .worldline(&ship.view(), "ship", Style::series(1))
        .unwrap();
    let all = elements(&diagram);
    let lines = find(&all, "worldline", "polyline");
    assert_eq!(lines.len(), 2);
    let layout = diagram.layout();

    // straight worldlines take two points each
    let rest_points = lines[0].points();
    assert_eq!(rest_points.len(), 2);
    assert_eq!(rest_points[0].0, rest_points[1].0, "vertical");
    assert!((rest_points[0].0 - layout.to_pixels(1.5, 0.0).0).abs() <= ROUNDING);
    // it runs the whole height of the plot and no further
    assert!((rest_points[0].1 - (layout.plot_top + layout.plot_height)).abs() <= ROUNDING);
    assert!((rest_points[1].1 - layout.plot_top).abs() <= ROUNDING);

    let ship_points = lines[1].points();
    assert_eq!(ship_points.len(), 2);
    let (dx, dy) = (
        ship_points[1].0 - ship_points[0].0,
        ship_points[1].1 - ship_points[0].1,
    );
    // time over distance is 1 / v; pixel y grows downwards
    assert!((-dy / dx - 1.0 / 0.6).abs() < 1e-3, "slope {}", -dy / dx);
    // and it passes through the origin
    let through = layout.to_pixels(0.0, 0.0);
    assert!(distance_to_polylines(through, &[ship_points]) <= ROUNDING);
}

#[test]
fn boosting_into_an_observers_rest_frame_stands_their_worldline_up() {
    let ship = moving(0.6, -100.0, 100.0);
    let station = at_rest(0.0, -100.0, 100.0);
    let frame = Frame::moving_with(&ship.view(), 3.0, Vec3::X).unwrap();
    assert!((frame.rapidity - 0.6_f64.atanh()).abs() < 1e-15);
    assert!(frame.origin.t == 3.0 && frame.origin.x.dist(Vec3::new(1.8, 0.0, 0.0)) < 1e-12);

    let mut diagram = Diagram::new(square_view().frame(frame)).unwrap();
    diagram
        .worldline(&ship.view(), "ship", Style::series(0))
        .unwrap()
        .worldline(&station.view(), "station", Style::series(1))
        .unwrap()
        .light_cone(on_axis(3.0, 1.8), Cone::Both, "light", Style::light())
        .unwrap();
    let all = elements(&diagram);
    let lines = find(&all, "worldline", "polyline");

    // the ship is at rest in its own frame, at the position of the frame's fixed event
    let ship_points = lines[0].points();
    assert_eq!(ship_points.len(), 2);
    assert!((ship_points[0].0 - ship_points[1].0).abs() <= ROUNDING);
    let layout = diagram.layout();
    assert!((ship_points[0].0 - layout.to_pixels(1.8, 0.0).0).abs() <= ROUNDING);

    // the station now moves the other way at the same speed
    let station_points = lines[1].points();
    let (dx, dy) = (
        station_points[1].0 - station_points[0].0,
        station_points[1].1 - station_points[0].1,
    );
    assert!((-dy / dx + 1.0 / 0.6).abs() < 1e-3, "slope {}", -dy / dx);

    // and light is still at 45 degrees
    let cone = find(&all, "cone", "line");
    assert_eq!(cone.len(), 4);
    for line in cone {
        let (dx, dy) = run(line);
        assert!((dx.abs() - dy.abs()).abs() <= 2.0 * ROUNDING);
    }

    // an accelerating ship is momentarily upright in the frame that moves with it
    let rocket = single(Segment::hyperbolic(
        0.0,
        50.0,
        Vec3::ZERO,
        Vec3::X,
        0.4,
        0.0,
    ));
    let rocket = rocket.view();
    let comoving = Frame::moving_with(&rocket, 2.0, Vec3::X).unwrap();
    let diagram =
        Diagram::new(View::equal_scales((-4.0, 6.0), (-2.0, 6.0), 60.0).frame(comoving)).unwrap();
    let before = diagram.coordinates(Event::new(1.99, rocket.at(1.99).x));
    let after = diagram.coordinates(Event::new(2.01, rocket.at(2.01).x));
    let here = diagram.coordinates(Event::new(2.0, rocket.at(2.0).x));
    // on both sides of the moment the ship is ahead of where it rests: the worldline touches the vertical
    assert!(before.0 > here.0 && after.0 > here.0);
    assert!(before.0 - here.0 < 1e-4 && after.0 - here.0 < 1e-4);
}

#[test]
fn a_boost_of_zero_rapidity_changes_nothing() {
    let ship = moving(0.3, -100.0, 100.0);
    let draw = |view: View| {
        let mut diagram = Diagram::new(view).unwrap();
        diagram
            .worldline(&ship.view(), "ship", Style::series(0))
            .unwrap()
            .light_cone(on_axis(1.0, 1.0), Cone::Both, "light", Style::light())
            .unwrap();
        diagram.render()
    };
    let lab = draw(square_view());
    let still = draw(square_view().frame(Frame::new(on_axis(7.0, -3.0), 0.0)));
    assert_eq!(lab, still);
}

#[test]
fn proper_time_ticks_are_gamma_dtau_apart_in_coordinate_time() {
    // gamma = 5/3 at 0.8 c
    let clock = moving(0.8, -100.0, 100.0);
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .proper_time_ticks(
            &clock.view(),
            0.75,
            TickNumbers::Right,
            "ticks",
            Style::series(0),
        )
        .unwrap();
    let all = elements(&diagram);
    let ticks = find(&all, "ticks", "circle");
    // the view is 12 units of time high: ticks at tau = -2.25 ... 4.5 fit, 1.25 apart in t
    assert_eq!(ticks.len(), 10);
    let layout = diagram.layout();
    let scale = layout.plot_height / (layout.t.1 - layout.t.0);
    for pair in ticks.windows(2) {
        let (a, b) = (pair[0].centre(), pair[1].centre());
        // later ticks are higher on the page
        assert!(((a.1 - b.1) - 0.75 * (5.0 / 3.0) * scale).abs() <= 2.0 * ROUNDING);
        assert!(((b.0 - a.0) - 0.8 * 0.75 * (5.0 / 3.0) * scale).abs() <= 2.0 * ROUNDING);
    }
    // one of them is the origin, where the clock reads zero
    assert!(ticks
        .iter()
        .any(|t| close(t.centre(), layout.to_pixels(0.0, 0.0), ROUNDING)));
    // and the numbers beside them are the clock readings; the last tick is on the right-hand edge of the
    // plot, where its reading would run off it, so that one is left out
    let numbers: Vec<&str> = find(&all, "labels", "text")
        .iter()
        .filter(|t| t.attr("stroke").is_none())
        .map(|t| t.text.as_str())
        .collect();
    assert_eq!(
        numbers,
        [
            "\u{2212}2.25",
            "\u{2212}1.50",
            "\u{2212}0.75",
            "0.00",
            "0.75",
            "1.50",
            "2.25",
            "3.00",
            "3.75"
        ]
    );
    assert!(close(
        ticks[9].centre(),
        (
            layout.plot_left + layout.plot_width,
            layout.to_pixels(0.0, 7.5).1
        ),
        ROUNDING
    ));

    // in the clock's own frame the same ticks are dtau apart and stacked vertically
    let own = Frame::moving_with(&clock.view(), 0.0, Vec3::X).unwrap();
    let mut diagram = Diagram::new(square_view().frame(own)).unwrap();
    diagram
        .proper_time_ticks(
            &clock.view(),
            0.75,
            TickNumbers::None,
            "ticks",
            Style::series(0),
        )
        .unwrap();
    let all = elements(&diagram);
    let ticks = find(&all, "ticks", "circle");
    // tau = -3.75 ... 7.5
    assert_eq!(ticks.len(), 16);
    for pair in ticks.windows(2) {
        let (a, b) = (pair[0].centre(), pair[1].centre());
        assert!(((a.1 - b.1) - 0.75 * scale).abs() <= 2.0 * ROUNDING);
        assert!((a.0 - b.0).abs() <= 2.0 * ROUNDING);
    }
}

#[test]
fn ticks_on_an_accelerating_clock_follow_its_proper_time() {
    let rocket = single(Segment::hyperbolic(
        0.0,
        50.0,
        Vec3::ZERO,
        Vec3::X,
        0.5,
        0.0,
    ));
    let rocket = rocket.view();
    let mut diagram = Diagram::new(View::equal_scales((-1.0, 9.0), (0.0, 10.0), 50.0)).unwrap();
    diagram
        .proper_time_ticks(&rocket, 0.5, TickNumbers::None, "", Style::series(0))
        .unwrap();
    let all = elements(&diagram);
    let ticks = find(&all, "ticks", "circle");
    assert!(ticks.len() >= 8);
    for (k, tick) in ticks.iter().enumerate() {
        // x = (cosh(a tau) - 1) / a, t = sinh(a tau) / a
        let tau = 0.5 * k as f64;
        let expected = on_axis((0.5 * tau).sinh() / 0.5, ((0.5 * tau).cosh() - 1.0) / 0.5);
        assert!(
            close(tick.centre(), pixel(&diagram, expected), ROUNDING),
            "tick {k}"
        );
    }
    // the next tick would be off the top or the side
    let tau = 0.5 * ticks.len() as f64;
    assert!((0.5 * tau).sinh() / 0.5 > 10.0 || ((0.5 * tau).cosh() - 1.0) / 0.5 > 9.0);
}

fn boost_coast_brake() -> WorldlineBuf {
    // from rest at the origin: a unit of proper time boosting, a coast, a brake to rest, then stay
    let boost = Segment::hyperbolic(0.0, 1.0_f64.sinh(), Vec3::ZERO, Vec3::X, 1.0, 0.0);
    let s1 = boost.at_tau(1.0);
    let coast = Segment::inertial(s1.t, s1.t + 3.0, s1.x, s1.v);
    let s2 = coast.at(s1.t + 3.0);
    let brake = Segment::hyperbolic(
        s1.t + 3.0,
        s1.t + 3.0 + 1.0_f64.sinh(),
        s2.x,
        Vec3::X,
        -1.0,
        1.0,
    );
    let s3 = brake.at_tau(1.0);
    WorldlineBuf::new(
        [
            Segment::rest(-50.0, 0.0, Vec3::ZERO),
            boost,
            coast,
            brake,
            Segment::rest(s3.t, 50.0, s3.x),
        ],
        0.0,
        0.0,
    )
    .unwrap()
}

#[test]
fn what_is_seen_is_a_light_ray_that_ends_on_the_source() {
    let source = boost_coast_brake();
    let source = source.view();
    for (frame, observer) in [
        (None, on_axis(7.0, -1.0)),
        (None, on_axis(3.5, 5.5)),
        (Some(Frame::new(origin(), 0.5)), on_axis(5.0, -0.5)),
    ] {
        let mut view = View::equal_scales((-4.0, 6.0), (-2.0, 8.0), 60.0);
        view.frame = frame;
        let mut diagram = Diagram::new(view).unwrap();
        diagram
            .worldline(&source, "source", Style::series(0))
            .unwrap()
            .seen(observer, &source, "seen", Style::light())
            .unwrap();
        let all = elements(&diagram);
        let lines = find(&all, "seen", "line");
        assert_eq!(lines.len(), 1);
        let (from, to) = lines[0].ends();

        // the ray is null: 45 degrees on the page
        let (dx, dy) = run(lines[0]);
        assert!((dx.abs() - dy.abs()).abs() <= 2.0 * ROUNDING, "{dx} {dy}");
        assert!(dy < 0.0, "light travels up the page");

        // it arrives at the observer and left from the event Worldline::send finds
        assert!(close(to, pixel(&diagram, observer), ROUNDING));
        let emitted = source.send(observer);
        let emission = Event::new(emitted, source.at(emitted).x);
        assert!((observer.t - emitted - emission.x.dist(observer.x)).abs() < 1e-12);
        assert!(close(from, pixel(&diagram, emission), ROUNDING));

        // that event is on the worldline as drawn, and carries the ring
        let drawn: Vec<Vec<(f64, f64)>> = find(&all, "worldline", "polyline")
            .iter()
            .map(|p| p.points())
            .collect();
        assert!(distance_to_polylines(from, &drawn) <= CURVE_TOLERANCE + ROUNDING);
        let rings = find(&all, "seen-source", "circle");
        assert_eq!(rings.len(), 1);
        assert!(close(rings[0].centre(), from, ROUNDING));
    }
}

/// The largest distance, in pixels, from the true worldline to the polyline drawn for it, over `samples`
/// moments between `t0` and `t1`; the number of points drawn; and the number of samples that were in view.
///
/// A sample within 0.15 px of the edge of the plot is passed over. Where a curve leaves the plot the chord
/// that stands for it is cut at the edge, and the curve can run on for the width of the tolerance before it
/// crosses too; that sliver lies under the frame of the plot.
fn deviation(
    view: View,
    worldline: &Worldline<'_>,
    t0: f64,
    t1: f64,
    samples: u32,
) -> (f64, usize, usize) {
    let mut diagram = Diagram::new(view).unwrap();
    diagram
        .worldline(worldline, "curve", Style::series(0))
        .unwrap();
    let all = elements(&diagram);
    let drawn: Vec<Vec<(f64, f64)>> = find(&all, "worldline", "polyline")
        .iter()
        .map(|p| p.points())
        .collect();
    let layout = diagram.layout();
    let (right, bottom) = (
        layout.plot_left + layout.plot_width,
        layout.plot_top + layout.plot_height,
    );
    let mut worst: f64 = 0.0;
    let mut inside = 0;
    for i in 0..=samples {
        let t = t0 + (t1 - t0) * f64::from(i) / f64::from(samples);
        let (x, time) = diagram.coordinates(Event::new(t, worldline.at(t).x));
        let p = layout.to_pixels(x, time);
        let margin = (p.0 - layout.plot_left)
            .min(right - p.0)
            .min(p.1 - layout.plot_top)
            .min(bottom - p.1);
        if margin.is_nan() || margin < 0.15 {
            continue;
        }
        inside += 1;
        worst = worst.max(distance_to_polylines(p, &drawn));
    }
    (worst, drawn.iter().map(Vec::len).sum(), inside)
}

#[test]
fn curved_worldlines_stay_within_a_quarter_pixel() {
    let view = View::equal_scales((-1.0, 9.0), (0.0, 10.0), 64.0);
    let burn = {
        let g = Burn::new(0.0, f64::NAN, Vec3::ZERO, Vec3::X, 0.0, 0.5, 0.8).unwrap();
        let mut segment = Segment::Burn(g);
        segment.set_t1(segment.t_of_tau(0.9 * g.t_burn));
        segment
    };
    let steered = Segment::Steered(
        Steered::new(
            0.0,
            20.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.2, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
            0.7,
        )
        .unwrap(),
    );
    let cases: [(&str, Segment, View); 7] = [
        (
            "from rest at 0.8",
            Segment::hyperbolic(0.0, 20.0, Vec3::ZERO, Vec3::X, 0.8, 0.0),
            view,
        ),
        (
            "braking through a turnaround",
            Segment::hyperbolic(0.0, 20.0, Vec3::ZERO, Vec3::X, -0.5, 1.2),
            view,
        ),
        (
            "a hard burn that leaves the view",
            Segment::hyperbolic(0.0, 20.0, Vec3::ZERO, Vec3::X, 40.0, 0.0),
            view,
        ),
        (
            "a circle seen edge on",
            Segment::circular(0.0, 40.0, Vec3::ZERO, 3.0, 0.3, Vec3::X, Vec3::Y, 0.4),
            View::equal_scales((-4.0, 4.0), (0.0, 40.0), 24.0),
        ),
        (
            "a burn that spends mass",
            burn,
            View::equal_scales((-0.2, 1.8), (0.0, 2.0), 300.0),
        ),
        ("thrust at an angle", steered, view),
        (
            "from rest at 0.8, seen from a moving frame",
            Segment::hyperbolic(0.0, 20.0, Vec3::ZERO, Vec3::X, 0.8, 0.0),
            View::equal_scales((-6.0, 4.0), (0.0, 10.0), 64.0).frame(Frame::new(origin(), 0.8)),
        ),
    ];
    for (name, segment, view) in cases {
        let worldline = single(segment);
        let (worst, points, inside) =
            deviation(view, &worldline.view(), segment.t0(), segment.t1(), 6000);
        assert!(inside > 1000, "{name}: only {inside} samples in view");
        assert!(
            worst <= CURVE_TOLERANCE + 0.008,
            "{name}: {worst:.3} px off with {points} points"
        );
        assert!((3..600).contains(&points), "{name}: {points} points");
    }
}

/// A small generator of repeatable numbers in `0..1`, so that the sweep below needs no dependency.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn between(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
}

#[test]
fn curves_of_every_kind_stay_within_tolerance_in_any_view() {
    let mut numbers = Numbers(2026);
    let (mut worst_of_all, mut most_points, mut drawn): (f64, usize, usize) = (0.0, 0, 0);
    for case in 0..400 {
        let t0 = numbers.between(-5.0, 5.0);
        let p0 = Vec3::new(numbers.between(-3.0, 3.0), numbers.between(-1.0, 1.0), 0.0);
        let segment = match case % 5 {
            0 => Segment::hyperbolic(
                t0,
                t0 + 30.0,
                p0,
                Vec3::X,
                numbers.between(0.01, 5.0),
                numbers.between(-2.0, 2.0),
            ),
            // braking hard, through a turnaround
            1 => Segment::hyperbolic(
                t0,
                t0 + 30.0,
                p0,
                Vec3::X,
                -numbers.between(0.01, 50.0),
                numbers.between(0.0, 3.0),
            ),
            2 => {
                let (r, beta) = (numbers.between(0.05, 6.0), numbers.between(0.01, 0.99));
                let phase = numbers.between(0.0, 6.0);
                Segment::circular(t0, t0 + 60.0, p0, r, beta / r, Vec3::X, Vec3::Y, phase)
            }
            3 => {
                let thrust = Vec3::new(numbers.between(-1.0, 1.0), numbers.between(-1.0, 1.0), 0.1);
                let v = Vec3::new(numbers.between(-0.9, 0.9), numbers.between(-0.3, 0.3), 0.0);
                let alpha = numbers.between(0.05, 4.0);
                Segment::Steered(Steered::new(t0, t0 + 30.0, p0, v, thrust, alpha).unwrap())
            }
            _ => {
                let sign = if numbers.next() < 0.5 { -1.0 } else { 1.0 };
                let alpha = sign * numbers.between(0.05, 3.0);
                let (phi, w) = (numbers.between(-1.0, 1.0), numbers.between(0.05, 1.0));
                let g = Burn::new(t0, f64::NAN, p0, Vec3::X, phi, alpha, w).unwrap();
                let mut segment = Segment::Burn(g);
                segment.set_t1(segment.t_of_tau(0.95 * g.t_burn));
                segment
            }
        };
        // a window of any shape somewhere near the curve, half the time in a moving frame
        let (cx, ct) = (numbers.between(-4.0, 6.0), t0 + numbers.between(0.0, 12.0));
        let (half_x, half_t) = (numbers.between(0.5, 12.0), numbers.between(0.5, 12.0));
        let mut view = View::new(
            (cx - half_x, cx + half_x),
            (ct - half_t, ct + half_t),
            numbers.between(100.0, 1200.0),
            numbers.between(100.0, 1200.0),
        );
        if numbers.next() < 0.5 {
            let pivot = on_axis(numbers.between(-3.0, 3.0), numbers.between(-3.0, 3.0));
            view = view.frame(Frame::new(pivot, numbers.between(-1.5, 1.5)));
        }
        let worldline = single(segment);
        let (worst, points, inside) =
            deviation(view, &worldline.view(), segment.t0(), segment.t1(), 3000);
        if inside > 2 {
            drawn += 1;
            assert!(points >= 2, "case {case}: in view but not drawn");
            assert!(
                worst <= CURVE_TOLERANCE + 0.008,
                "case {case}: {worst:.3} px off with {points} points in {view:?}"
            );
            worst_of_all = worst_of_all.max(worst);
        }
        most_points = most_points.max(points);
    }
    assert!(drawn > 250, "only {drawn} curves crossed their view");
    // the sampler is not wasteful either
    assert!(most_points < 1500, "{most_points} points for one curve");
    assert!(
        worst_of_all > 0.01,
        "the comparison saw no curvature at all"
    );
}

#[test]
fn nothing_is_drawn_past_the_end_of_a_worldline_or_the_edge_of_the_plot() {
    // A trip that exists only from t = 1 to t = 5, in a view from -2 to 10.
    let trip = single(Segment::hyperbolic(1.0, 5.0, Vec3::ZERO, Vec3::X, 0.3, 0.0));
    let trip = trip.view();
    let view = View::equal_scales((-2.0, 6.0), (-2.0, 10.0), 50.0);
    let mut diagram = Diagram::new(view).unwrap();
    diagram
        .worldline(&trip, "trip", Style::series(0))
        .unwrap()
        .proper_time_ticks(&trip, 0.5, TickNumbers::None, "ticks", Style::series(0))
        .unwrap();
    let all = elements(&diagram);
    let lines = find(&all, "worldline", "polyline");
    assert_eq!(lines.len(), 1);
    let points = lines[0].points();
    assert!(close(
        points[0],
        pixel(&diagram, trip_event(&trip, 1.0)),
        ROUNDING
    ));
    assert!(close(
        *points.last().unwrap(),
        pixel(&diagram, trip_event(&trip, 5.0)),
        ROUNDING
    ));
    // ticks only between the two ends: tau runs from 0 to a little under 3.3
    let ticks = find(&all, "ticks", "circle");
    assert_eq!(ticks.len(), 1 + (trip.tau(5.0) / 0.5).floor() as usize);
    let (top, bottom) = (points.last().unwrap().1, points[0].1);
    assert!(ticks
        .iter()
        .all(|t| t.centre().1 >= top - ROUNDING && t.centre().1 <= bottom + ROUNDING));

    // only part of it, on request
    let mut diagram = Diagram::new(view).unwrap();
    diagram
        .worldline_between(&trip, 2.0, 3.0, "part", Style::series(0))
        .unwrap();
    let all = elements(&diagram);
    let points = find(&all, "worldline", "polyline")[0].points();
    assert!(close(
        points[0],
        pixel(&diagram, trip_event(&trip, 2.0)),
        ROUNDING
    ));
    assert!(close(
        *points.last().unwrap(),
        pixel(&diagram, trip_event(&trip, 3.0)),
        ROUNDING
    ));

    // A ship that flies out of the side of the view, with light and a grid: every mark stays in the plot.
    let ship = single(Segment::hyperbolic(
        0.0,
        50.0,
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::X,
        0.7,
        0.0,
    ));
    let frame = Frame::from_velocity(origin(), 0.4).unwrap();
    let mut diagram = Diagram::new(view).unwrap();
    diagram
        .worldline(&ship.view(), "ship", Style::series(0))
        .unwrap()
        .light_cone(on_axis(9.0, 5.0), Cone::Both, "light", Style::light())
        .unwrap()
        .frame_grid(frame, 1.0, "grid", Style::guide())
        .unwrap()
        .signal(
            on_axis(-5.0, -8.0),
            on_axis(15.0, 12.0),
            "ray",
            Style::light(),
        )
        .unwrap()
        .rindler_horizon(&ship.segments()[0].seg, "horizon", Style::horizon())
        .unwrap()
        // events outside the view are not drawn at all
        .event(on_axis(20.0, 0.0), "late", Style::series(1))
        .unwrap()
        .event(on_axis(0.0, 40.0), "far", Style::series(1))
        .unwrap();
    let layout = diagram.layout();
    let (left, top) = (layout.plot_left - ROUNDING, layout.plot_top - ROUNDING);
    let (right, bottom) = (
        layout.plot_left + layout.plot_width + ROUNDING,
        layout.plot_top + layout.plot_height + ROUNDING,
    );
    let inside = |p: (f64, f64)| p.0 >= left && p.0 <= right && p.1 >= top && p.1 <= bottom;
    let all = elements(&diagram);
    let mut checked = 0;
    for class in ["worldline", "cone", "frame-grid", "signal", "horizon"] {
        for line in find(&all, class, "line") {
            let (a, b) = line.ends();
            assert!(inside(a) && inside(b), "{class} line {a:?} {b:?}");
            checked += 1;
        }
        for line in find(&all, class, "polyline") {
            assert!(line.points().into_iter().all(inside), "{class} polyline");
            checked += 1;
        }
    }
    assert!(checked > 20);
    assert!(find(&all, "event", "circle").is_empty());
    // the ship's polyline ends exactly on the right-hand edge, where it leaves
    let ship_points = find(&all, "worldline", "polyline")[0].points();
    assert!((ship_points.last().unwrap().0 - (right - ROUNDING)).abs() <= ROUNDING);
}

fn trip_event(worldline: &Worldline<'_>, t: f64) -> Event {
    Event::new(t, worldline.at(t).x)
}

#[test]
fn another_frames_lines_tilt_by_its_speed() {
    let frame = Frame::from_velocity(origin(), 0.5).unwrap();
    let gamma = 1.0 / (1.0_f64 - 0.25).sqrt();
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .simultaneity(frame, 2.0, "now", Style::guide())
        .unwrap()
        .constant_position(frame, -1.0, "here", Style::guide())
        .unwrap()
        .frame_axes(frame, "axes", Style::guide())
        .unwrap();
    let all = elements(&diagram);
    let layout = diagram.layout();

    // a line of simultaneity rises by v for each unit of position
    let now = find(&all, "simultaneity", "line");
    assert_eq!(now.len(), 1);
    let (dx, dy) = run(now[0]);
    assert!((-dy / dx - 0.5).abs() < 1e-3);
    // it holds the event with t' = 2, x' = 0: the lab event (gamma 2, gamma v 2)
    let on_it = layout.to_pixels(gamma * 0.5 * 2.0, gamma * 2.0);
    assert!(
        distance_to_polylines(on_it, &[vec![now[0].ends().0, now[0].ends().1]]) <= 2.0 * ROUNDING
    );

    // a line of constant position leans by v for each unit of time
    let here = find(&all, "position", "line");
    assert_eq!(here.len(), 1);
    let (dx, dy) = run(here[0]);
    assert!((dx / -dy - 0.5).abs() < 1e-3);
    // it holds the event with t' = 0, x' = -1: the lab event (-gamma v, -gamma)
    let on_it = layout.to_pixels(-gamma, -gamma * 0.5);
    assert!(
        distance_to_polylines(on_it, &[vec![here[0].ends().0, here[0].ends().1]]) <= 2.0 * ROUNDING
    );

    // the axes pass through the frame's origin, one of each kind, mirrored in the light cone
    let axes = find(&all, "frame-axes", "line");
    assert_eq!(axes.len(), 2);
    let through = layout.to_pixels(0.0, 0.0);
    let slopes: Vec<f64> = axes
        .iter()
        .map(|line| {
            let (a, b) = line.ends();
            assert!(distance_to_polylines(through, &[vec![a, b]]) <= 2.0 * ROUNDING);
            let (dx, dy) = run(line);
            -dy / dx
        })
        .collect();
    assert!((slopes[0] - 0.5).abs() < 1e-3 && (slopes[1] - 2.0).abs() < 1e-3);

    // Drawn in that frame its own grid is square, and the lab's grid tilts the other way.
    let mut diagram = Diagram::new(square_view().frame(frame)).unwrap();
    diagram
        .frame_grid(frame, 2.5, "own", Style::guide())
        .unwrap()
        .simultaneity(Frame::new(origin(), 0.0), 0.0, "lab now", Style::guide())
        .unwrap();
    let all = elements(&diagram);
    let grid = find(&all, "frame-grid", "line");
    // t' = -2.5 ... 7.5 and x' = -5 ... 5, every 2.5
    assert_eq!(grid.len(), 10);
    assert!(grid.iter().all(|line| {
        let (dx, dy) = run(line);
        dx == 0.0 || dy == 0.0
    }));
    let (dx, dy) = run(find(&all, "simultaneity", "line")[0]);
    assert!((-dy / dx + 0.5).abs() < 1e-3);
}

#[test]
fn the_rindler_horizon_is_the_light_ray_through_the_pivot() {
    let drive = Hyperbolic::new(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, 0.5, 0.0);
    let horizon = RindlerHorizon::of(&drive);
    assert_eq!(horizon.pivot, on_axis(0.0, -2.0));
    for frame in [None, Some(Frame::new(origin(), -0.6))] {
        let mut view = square_view();
        view.frame = frame;
        let mut diagram = Diagram::new(view).unwrap();
        diagram
            .rindler_horizon(&Segment::Hyperbolic(drive), "horizon", Style::horizon())
            .unwrap();
        let all = elements(&diagram);
        let lines = find(&all, "horizon", "line");
        assert_eq!(lines.len(), 1);
        let (a, b) = lines[0].ends();
        // null, rising to the right, through the pivot
        let (dx, dy) = run(lines[0]);
        assert!((dx + dy).abs() <= 2.0 * ROUNDING && dx > 0.0);
        assert!(
            distance_to_polylines(pixel(&diagram, horizon.pivot), &[vec![a, b]]) <= 2.0 * ROUNDING
        );
    }

    // In the lab frame: the events on the drawn line are exactly those whose light just fails to arrive.
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .rindler_horizon(&Segment::Hyperbolic(drive), "horizon", Style::horizon())
        .unwrap();
    let all = elements(&diagram);
    let layout = diagram.layout();
    let (a, b) = find(&all, "horizon", "line")[0].ends();
    for p in [a, b] {
        let (x, t) = coordinates(&layout, p);
        assert!(horizon.margin(on_axis(t, x)).abs() < 1e-3);
    }

    // A ship pushed backwards has its horizon ahead of it, falling to the right.
    let reverse = Segment::hyperbolic(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, -0.5, 0.0);
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .rindler_horizon(&reverse, "horizon", Style::horizon())
        .unwrap();
    let all = elements(&diagram);
    let (dx, dy) = run(find(&all, "horizon", "line")[0]);
    assert!((dx - dy).abs() <= 2.0 * ROUNDING);

    // Thrust at an angle fixed aboard has the same horizon when the angle is zero.
    let steered = Segment::Steered(
        Steered::new(0.0, f64::INFINITY, Vec3::ZERO, Vec3::ZERO, Vec3::X, 0.5).unwrap(),
    );
    let draw = |segment: &Segment| {
        let mut diagram = Diagram::new(square_view()).unwrap();
        diagram
            .rindler_horizon(segment, "horizon", Style::horizon())
            .unwrap();
        diagram.render()
    };
    assert_eq!(draw(&steered), draw(&Segment::Hyperbolic(drive)));

    // Anything else has none.
    let mut diagram = Diagram::new(square_view()).unwrap();
    for segment in [
        Segment::rest(0.0, 1.0, Vec3::ZERO),
        Segment::inertial(0.0, 1.0, Vec3::ZERO, Vec3::new(0.5, 0.0, 0.0)),
        Segment::null(0.0, 1.0, Vec3::ZERO, Vec3::X),
    ] {
        assert_eq!(
            diagram
                .rindler_horizon(&segment, "", Style::horizon())
                .err(),
            Some(Error::NoHorizon)
        );
    }
}

#[test]
fn signals_must_be_light_like() {
    let mut diagram = Diagram::new(square_view()).unwrap();
    let style = Style::light();
    // light, along the axis and at an angle to it
    assert!(diagram
        .signal(origin(), on_axis(3.0, 3.0), "", style)
        .is_ok());
    assert!(diagram
        .signal(origin(), on_axis(3.0, -3.0), "", style)
        .is_ok());
    assert!(diagram
        .signal(
            origin(),
            Event::new(5.0, Vec3::new(3.0, 4.0, 0.0)),
            "",
            style
        )
        .is_ok());
    // slower than light, faster than light, and backwards in time
    for to in [
        on_axis(3.0, 2.9),
        on_axis(3.0, 3.1),
        on_axis(-3.0, 3.0),
        on_axis(3.0, 0.0),
    ] {
        assert_eq!(
            diagram.signal(origin(), to, "", style).err(),
            Some(Error::NotLightLike),
            "{to:?}"
        );
    }
    // an arrival found by the core's light-cone solver passes
    let ship = single(Segment::hyperbolic(
        0.0,
        50.0,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::X,
        0.2,
        0.0,
    ));
    let from = on_axis(0.5, 0.0);
    let arrival = ship.view().recv(from);
    let to = Event::new(arrival, ship.view().at(arrival).x);
    assert!(diagram.signal(from, to, "", style).is_ok());

    // The arrowhead sits half way and points the way the light goes.
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .signal(origin(), on_axis(4.0, -4.0), "", style)
        .unwrap();
    let all = elements(&diagram);
    let line = find(&all, "signal", "line")[0].ends();
    let arrow = find(&all, "signal", "polygon")[0].points();
    assert_eq!(arrow.len(), 3);
    let middle = (0.5 * (line.0 .0 + line.1 .0), 0.5 * (line.0 .1 + line.1 .1));
    let centre = (
        (arrow[0].0 + arrow[1].0 + arrow[2].0) / 3.0,
        (arrow[0].1 + arrow[1].1 + arrow[2].1) / 3.0,
    );
    assert!((centre.0 - middle.0).abs() < 2.0 && (centre.1 - middle.1).abs() < 2.0);
    // the tip is up and to the left of the middle, as the light is
    assert!(arrow[0].0 < middle.0 && arrow[0].1 < middle.1);
}

#[test]
fn bad_input_is_refused_without_panicking() {
    let view = square_view();
    let ship = moving(0.5, -100.0, 100.0);
    let ship = ship.view();
    let style = Style::series(0);
    let nan = f64::NAN;
    let mut diagram = Diagram::new(view).unwrap();

    assert_eq!(
        Diagram::new(View::new((1.0, 1.0), (0.0, 1.0), 100.0, 100.0)).err(),
        Some(Error::EmptyView)
    );
    assert_eq!(
        Diagram::new(View::new((0.0, 1.0), (0.0, nan), 100.0, 100.0)).err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        Diagram::new(view.axis(Vec3::ZERO, Vec3::ZERO)).err(),
        Some(Error::ZeroDirection)
    );
    assert_eq!(
        Diagram::new(View::new((0.0, 1.0), (0.0, 1.0), 1e9, 100.0)).err(),
        Some(Error::TooLarge)
    );

    for width in [0.0, -1.0, nan, f64::INFINITY, 65.0] {
        let bad = style.with_width(width);
        assert_eq!(
            diagram.worldline(&ship, "", bad).err(),
            Some(Error::InvalidStyle)
        );
        assert_eq!(
            diagram.event(origin(), "", bad).err(),
            Some(Error::InvalidStyle)
        );
        assert_eq!(
            diagram.light_cone(origin(), Cone::Both, "", bad).err(),
            Some(Error::InvalidStyle)
        );
    }
    let nowhere = Event::new(nan, Vec3::ZERO);
    let far = Event::new(0.0, Vec3::new(f64::INFINITY, 0.0, 0.0));
    for event in [nowhere, far] {
        assert_eq!(
            diagram.event(event, "", style).err(),
            Some(Error::NonFinite)
        );
        assert_eq!(
            diagram.light_cone(event, Cone::Past, "", style).err(),
            Some(Error::NonFinite)
        );
        assert_eq!(
            diagram.signal(event, origin(), "", style).err(),
            Some(Error::NonFinite)
        );
        assert_eq!(
            diagram.seen(event, &ship, "", style).err(),
            Some(Error::NonFinite)
        );
    }
    assert_eq!(
        diagram.worldline_between(&ship, nan, 1.0, "", style).err(),
        Some(Error::NonFinite)
    );

    // steps
    let frame = Frame::new(origin(), 0.3);
    for step in [0.0, -1.0] {
        assert_eq!(
            diagram
                .proper_time_ticks(&ship, step, TickNumbers::None, "", style)
                .err(),
            Some(Error::OutOfRange)
        );
        assert_eq!(
            diagram.frame_grid(frame, step, "", style).err(),
            Some(Error::OutOfRange)
        );
    }
    assert_eq!(
        diagram
            .proper_time_ticks(&ship, nan, TickNumbers::None, "", style)
            .err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        diagram.frame_grid(frame, f64::INFINITY, "", style).err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        diagram
            .proper_time_ticks(&ship, 1e-6, TickNumbers::None, "", style)
            .err(),
        Some(Error::TooManyMarks)
    );
    assert_eq!(
        diagram.frame_grid(frame, 1e-6, "", style).err(),
        Some(Error::TooManyMarks)
    );

    // frames
    let broken = Frame::new(origin(), nan);
    assert_eq!(
        diagram.simultaneity(broken, 0.0, "", style).err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        diagram.constant_position(frame, nan, "", style).err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        diagram.frame_axes(broken, "", style).err(),
        Some(Error::NonFinite)
    );
    assert_eq!(
        Frame::from_velocity(origin(), 1.5).err(),
        Some(Error::OutOfRange)
    );
    assert_eq!(
        Frame::moving_with(&ship, 0.0, Vec3::ZERO).err(),
        Some(Error::ZeroDirection)
    );
    assert_eq!(
        Frame::moving_with(&ship, nan, Vec3::X).err(),
        Some(Error::NonFinite)
    );

    // a circle with millions of turns in view cannot be drawn to tolerance
    let spinner = single(Segment::circular(
        -50.0,
        50.0,
        Vec3::ZERO,
        1e-7,
        5e6,
        Vec3::X,
        Vec3::Y,
        0.0,
    ));
    assert_eq!(
        diagram.worldline(&spinner.view(), "", style).err(),
        Some(Error::TooDetailed)
    );
    // a segment with a broken parameter gives positions that are not numbers
    let broken = single(Segment::hyperbolic(
        0.0,
        10.0,
        Vec3::ZERO,
        Vec3::X,
        nan,
        0.0,
    ));
    assert_eq!(
        diagram.worldline(&broken.view(), "", style).err(),
        Some(Error::NonFinite)
    );
    let coasting = Segment::hyperbolic(0.0, 10.0, Vec3::ZERO, Vec3::X, 0.0, 0.0);
    assert_eq!(
        diagram.rindler_horizon(&coasting, "", style).err(),
        Some(Error::NonFinite)
    );

    // nothing to see: the source did not exist yet, or sits behind the observer's past light cone for ever
    let late = at_rest(2.0, 0.0, 50.0);
    assert_eq!(
        diagram
            .seen(on_axis(1.0, 0.0), &late.view(), "", style)
            .err(),
        Some(Error::NothingSeen)
    );
    assert!(diagram
        .seen(on_axis(3.0, 0.0), &late.view(), "", style)
        .is_ok());
    let runaway = single(Segment::hyperbolic(
        0.0,
        f64::INFINITY,
        Vec3::ZERO,
        Vec3::X,
        1.0,
        0.0,
    ));
    assert_eq!(
        diagram
            .seen(on_axis(-3.0, 0.0), &runaway.view(), "", style)
            .err(),
        Some(Error::NothingSeen)
    );

    // None of the refused calls left anything behind.
    let empty = Diagram::new(view).unwrap();
    let mut expected = empty.clone();
    expected
        .seen(on_axis(3.0, 0.0), &late.view(), "", style)
        .unwrap();
    assert_eq!(diagram, expected);
}

#[test]
fn a_view_along_another_axis_projects_onto_it() {
    // A ship flying along y, drawn along y: the same picture as the ship along x drawn along x.
    let along_y = WorldlineBuf::new(
        [Segment::hyperbolic(
            0.0,
            50.0,
            Vec3::new(7.0, 1.0, 0.0),
            Vec3::Y,
            0.4,
            0.0,
        )],
        0.0,
        0.0,
    )
    .unwrap();
    let along_x = single(Segment::hyperbolic(
        0.0,
        50.0,
        Vec3::ZERO,
        Vec3::X,
        0.4,
        0.0,
    ));
    let draw = |view: View, worldline: &WorldlineBuf| {
        let mut diagram = Diagram::new(view).unwrap();
        diagram
            .worldline(&worldline.view(), "ship", Style::series(0))
            .unwrap()
            .proper_time_ticks(
                &worldline.view(),
                1.0,
                TickNumbers::Right,
                "",
                Style::series(0),
            )
            .unwrap();
        diagram.render()
    };
    let view = View::equal_scales((-1.0, 9.0), (0.0, 10.0), 40.0);
    // the direction need not be a unit vector
    let turned = view.axis(Vec3::new(7.0, 1.0, 0.0), Vec3::new(0.0, 3.0, 0.0));
    assert_eq!(draw(turned, &along_y), draw(view, &along_x));

    // A body circling in the xy plane is a sine wave along x: it swings between -r and +r.
    let circle = single(Segment::circular(
        0.0,
        100.0,
        Vec3::ZERO,
        2.0,
        0.25,
        Vec3::X,
        Vec3::Y,
        0.0,
    ));
    let mut diagram = Diagram::new(View::equal_scales((-3.0, 3.0), (0.0, 30.0), 20.0)).unwrap();
    diagram
        .worldline(&circle.view(), "orbit", Style::series(0))
        .unwrap();
    let all = elements(&diagram);
    let points = find(&all, "worldline", "polyline")[0].points();
    let layout = diagram.layout();
    let (min, max) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.0), hi.max(p.0))
        });
    assert!((min - layout.to_pixels(-2.0, 0.0).0).abs() < 0.3);
    assert!((max - layout.to_pixels(2.0, 0.0).0).abs() < 0.3);
}

#[test]
fn a_light_like_jump_is_a_straight_line_at_45_degrees() {
    // rest, a jump at the speed of light to x = 3, rest again
    let target = Vec3::new(3.0, 0.0, 0.0);
    let jump = Segment::Null(srkit_core::segment::Null::between(1.0, Vec3::ZERO, target));
    let traveller = WorldlineBuf::new(
        [
            Segment::rest(-1.0, 1.0, Vec3::ZERO),
            jump,
            Segment::rest(4.0, 7.0, target),
        ],
        0.0,
        0.0,
    )
    .unwrap();
    let mut diagram = Diagram::new(square_view()).unwrap();
    diagram
        .worldline(&traveller.view(), "traveller", Style::series(0))
        .unwrap()
        .proper_time_ticks(
            &traveller.view(),
            1.0,
            TickNumbers::None,
            "",
            Style::series(0),
        )
        .unwrap();
    let all = elements(&diagram);
    let points = find(&all, "worldline", "polyline")[0].points();
    // four corners: two points for each straight piece, shared at the joins
    assert_eq!(points.len(), 4);
    let layout = diagram.layout();
    for (point, (x, t)) in points
        .iter()
        .zip([(0.0, -1.0), (0.0, 1.0), (3.0, 4.0), (3.0, 7.0)])
    {
        assert!(close(*point, layout.to_pixels(x, t), ROUNDING));
    }
    // no proper time passes on the jump: the clock reads 1 on leaving and on arriving
    let ticks: Vec<(f64, f64)> = find(&all, "ticks", "circle")
        .iter()
        .map(|c| c.centre())
        .collect();
    let expected = [
        (0.0, -1.0),
        (0.0, 0.0),
        (3.0, 4.0),
        (3.0, 5.0),
        (3.0, 6.0),
        (3.0, 7.0),
    ];
    assert_eq!(ticks.len(), expected.len());
    for (tick, (x, t)) in ticks.iter().zip(expected) {
        assert!(
            close(*tick, layout.to_pixels(x, t), ROUNDING),
            "{tick:?} for ({x}, {t})"
        );
    }
}
