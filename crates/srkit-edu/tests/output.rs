//! Output: the SVG is well formed, self-contained, escaped, finite and the same bytes every time.
mod common;

use common::{at_rest, find, moving, parse, single, Element};
use srkit_core::{kinematics::Event, segment::Segment, trip::FlipBurn, Vec3};
use srkit_edu::{
    rindler_chase, rindler_chase_diagram, twin_paradox, twin_paradox_diagram, Colour, Cone, Dash,
    Diagram, Error, Frame, Placement, Style, Theme, TickNumbers, View,
};

fn on_axis(t: f64, x: f64) -> Event {
    Event::new(t, Vec3::new(x, 0.0, 0.0))
}

/// A diagram that uses every kind of mark, in a boosted frame.
fn everything(theme: Theme) -> Diagram {
    let origin = on_axis(0.0, 0.0);
    let station = at_rest(0.0, -50.0, 50.0);
    let ship = moving(0.6, -50.0, 50.0);
    let rocket = Segment::hyperbolic(
        0.0,
        f64::INFINITY,
        Vec3::new(-3.0, 0.0, 0.0),
        Vec3::X,
        0.5,
        0.0,
    );
    let orbit = single(Segment::circular(
        -50.0,
        50.0,
        Vec3::ZERO,
        4.0,
        0.2,
        Vec3::X,
        Vec3::Y,
        1.0,
    ));
    let frame = Frame::new(origin, 0.3);
    let mut diagram =
        Diagram::new(View::equal_scales((-6.0, 6.0), (-2.0, 8.0), 48.0).frame(frame)).unwrap();
    diagram
        .theme(theme)
        .title("Everything at once")
        .subtitle("Drawn in a frame moving at tanh 0.3")
        .axis_titles("x (ly)", "t (yr)");
    let observer = on_axis(6.0, 2.0);
    diagram
        .frame_grid(Frame::new(origin, 0.0), 2.0, "Lab grid", Style::guide())
        .unwrap()
        .frame_axes(
            Frame::new(origin, 0.0),
            "Lab axes",
            Style::guide().with_dash(Dash::Dashed),
        )
        .unwrap()
        .simultaneity(
            Frame::new(origin, 0.8),
            1.0,
            "A fast observer's now",
            Style::guide().with_dash(Dash::DashDotDot),
        )
        .unwrap()
        .constant_position(
            Frame::new(origin, 0.8),
            1.0,
            "A fast observer's here",
            Style::guide().with_dash(Dash::LongDash),
        )
        .unwrap()
        .light_cone(origin, Cone::Both, "Light cone", Style::light())
        .unwrap()
        .rindler_horizon(&rocket, "Horizon", Style::horizon())
        .unwrap()
        .worldline(&station.view(), "Station", Style::series(0))
        .unwrap()
        .worldline(&ship.view(), "Ship", Style::series(1))
        .unwrap()
        .worldline(&single(rocket).view(), "Rocket", Style::series(2))
        .unwrap()
        .worldline(&orbit.view(), "Orbit", Style::series(3))
        .unwrap()
        .proper_time_ticks(
            &ship.view(),
            1.0,
            TickNumbers::Right,
            "Ship clock",
            Style::series(1),
        )
        .unwrap()
        .proper_time_ticks(
            &station.view(),
            0.5,
            TickNumbers::Left,
            "Station clock",
            Style::series(0),
        )
        .unwrap()
        .signal(
            origin,
            on_axis(4.0, -4.0),
            "Signal",
            Style::light().with_colour(Colour::Rgb(200, 30, 30)),
        )
        .unwrap()
        .seen(
            observer,
            &station.view(),
            "What the observer sees",
            Style::light().with_width(2.0),
        )
        .unwrap()
        .event(origin, "Origin", Style::horizon())
        .unwrap()
        .event_placed(observer, "Observer", Placement::Above, Style::horizon())
        .unwrap();
    diagram
}

fn samples() -> Vec<(&'static str, String)> {
    let mut dark_twins = twin_paradox_diagram(4.2465, 1.0323).unwrap();
    let mut dark_chase = rindler_chase_diagram(1.0323).unwrap();
    let bare = Diagram::new(View::new((0.0, 1.0), (0.0, 1.0), 120.0, 80.0)).unwrap();
    vec![
        ("twins", twin_paradox(4.2465, 1.0323).unwrap()),
        ("twins, dark", dark_twins.theme(Theme::Dark).render()),
        ("twins, far and fast", twin_paradox(1e4, 30.0).unwrap()),
        ("twins, near and slow", twin_paradox(1e-3, 1e-3).unwrap()),
        ("chase", rindler_chase(1.0323).unwrap()),
        ("chase, dark", dark_chase.theme(Theme::Dark).render()),
        ("chase, gentle", rindler_chase(1e-6).unwrap()),
        ("everything", everything(Theme::Light).render()),
        ("everything, dark", everything(Theme::Dark).render()),
        ("nothing", bare.render()),
    ]
}

/// Attributes that hold one number.
const NUMERIC: [&str; 13] = [
    "x",
    "y",
    "x1",
    "y1",
    "x2",
    "y2",
    "cx",
    "cy",
    "r",
    "width",
    "height",
    "stroke-width",
    "font-size",
];

fn two_decimals(text: &str) -> bool {
    let digits = text.strip_prefix('-').unwrap_or(text);
    match digits.split_once('.') {
        Some((whole, fraction)) => {
            !whole.is_empty()
                && whole.bytes().all(|b| b.is_ascii_digit())
                && fraction.len() == 2
                && fraction.bytes().all(|b| b.is_ascii_digit())
        }
        None => false,
    }
}

#[test]
fn the_svg_is_well_formed_self_contained_and_finite() {
    for (name, svg) in samples() {
        let elements = parse(&svg).unwrap_or_else(|e| panic!("{name}: {e}"));
        let root = &elements[0];
        assert_eq!(root.name, "svg", "{name}");
        assert!(
            svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""),
            "{name}"
        );
        assert!(svg.ends_with("</svg>\n"), "{name}");
        assert_eq!(root.attr("xmlns"), Some("http://www.w3.org/2000/svg"));
        assert_eq!(
            root.attr("viewBox").unwrap(),
            format!(
                "0 0 {} {}",
                root.attr("width").unwrap(),
                root.attr("height").unwrap()
            )
        );
        assert_eq!(root.attr("role"), Some("img"));
        assert!(!root.attr("aria-label").unwrap().is_empty());
        assert_eq!(
            elements[1].name, "title",
            "{name}: the accessible name comes first"
        );
        assert_eq!(elements[2].name, "desc");

        // nothing that runs, and nothing fetched from elsewhere
        let allowed = [
            "svg", "title", "desc", "rect", "g", "line", "polyline", "polygon", "circle", "text",
            "tspan",
        ];
        for element in &elements {
            assert!(
                allowed.contains(&element.name.as_str()),
                "{name}: <{}>",
                element.name
            );
            for (attr, value) in &element.attrs {
                assert!(
                    !attr.starts_with("on") && !attr.contains("href") && attr != "style",
                    "{name}: {attr}"
                );
                assert!(
                    !value.contains("url(") && !value.contains("javascript:"),
                    "{name}: {attr}={value}"
                );
                if attr != "xmlns" {
                    assert!(!value.contains("://"), "{name}: {attr}={value}");
                }
            }
        }

        // every number is finite and has two decimals
        assert!(!svg.contains("NaN") && !svg.contains("inf"), "{name}");
        let (width, height) = (root.num("width"), root.num("height"));
        for element in &elements {
            for attr in NUMERIC {
                if let Some(value) = element.attr(attr) {
                    assert!(
                        two_decimals(value),
                        "{name}: <{}> {attr}={value}",
                        element.name
                    );
                    assert!(value.parse::<f64>().unwrap().is_finite());
                }
            }
            if element.attr("points").is_some() {
                for pair in element.attr("points").unwrap().split(' ') {
                    let (x, y) = pair.split_once(',').unwrap();
                    assert!(two_decimals(x) && two_decimals(y), "{name}: point {pair}");
                }
                for (x, y) in element.points() {
                    assert!(
                        (0.0..=width).contains(&x) && (0.0..=height).contains(&y),
                        "{name}: ({x}, {y})"
                    );
                }
            }
        }
        // the image is a whole number of pixels
        assert_eq!(width, width.round());
        assert_eq!(height, height.round());
    }
}

#[test]
fn labels_are_escaped() {
    let nasty = "a<b & \"c\" > 'd' </text><script>alert(1)</script>";
    let ship = moving(0.5, -50.0, 50.0);
    let mut diagram = Diagram::new(View::equal_scales((-3.0, 3.0), (0.0, 6.0), 50.0)).unwrap();
    diagram
        .title(nasty)
        .subtitle(nasty)
        .axis_titles(nasty, nasty);
    diagram
        .worldline(&ship.view(), nasty, Style::series(0))
        .unwrap()
        .event(on_axis(1.0, 0.5), nasty, Style::horizon())
        .unwrap()
        .proper_time_ticks(
            &ship.view(),
            1.0,
            TickNumbers::Right,
            "tab\tand\u{1}control",
            Style::series(0),
        )
        .unwrap();
    let svg = diagram.render();
    assert!(!svg.contains("<script") && !svg.contains("\"c\""));
    assert!(svg.contains("a&lt;b &amp; &quot;c&quot; &gt; &apos;d&apos; &lt;/text&gt;"));
    assert!(!svg.contains('\u{1}') && !svg.contains('\t'));

    // read back, every label is the text that went in
    let elements = parse(&svg).unwrap();
    assert_eq!(elements[0].attr("aria-label"), Some(nasty));
    let text_of = |class: &str| -> Vec<&Element> { find(&elements, class, "text") };
    assert_eq!(text_of("title")[0].text, nasty);
    assert_eq!(text_of("subtitle")[0].text, nasty);
    assert_eq!(text_of("axis-title").len(), 2);
    assert!(text_of("axis-title").iter().all(|t| t.text == nasty));
    assert!(text_of("legend").iter().any(|t| t.text == nasty));
    assert!(text_of("legend")
        .iter()
        .any(|t| t.text == "tab and\u{fffd}control"));
    assert!(text_of("labels").iter().any(|t| t.text == nasty));
    // the hover titles too
    assert!(elements
        .iter()
        .any(|e| e.name == "title" && e.ancestors.len() == 2 && e.text == nasty));
}

#[test]
fn the_same_calls_give_the_same_bytes() {
    let first = samples();
    let second = samples();
    assert_eq!(first.len(), second.len());
    for ((name, a), (_, b)) in first.iter().zip(&second) {
        assert!(a == b, "{name} differs between two runs");
    }
    // a diagram can be rendered twice, and a clone renders the same
    let diagram = everything(Theme::Light);
    assert_eq!(diagram.render(), diagram.render());
    assert_eq!(diagram.clone().render(), diagram.render());
}

/// The checked-in rendering of `twin_paradox(4.2465, 1.0323)`: Proxima Centauri and back at one gravity.
///
/// A change here is a change to what every caller gets. If it is meant, look at the new picture and write it
/// over the file: `cargo run -p srkit-edu --example twin_paradox > crates/srkit-edu/tests/golden/twin_paradox.svg`.
#[test]
fn the_twin_paradox_matches_its_golden_file() {
    let golden = include_str!("golden/twin_paradox.svg");
    let svg = twin_paradox(4.2465, 1.0323).unwrap();
    if svg != golden {
        let line = svg
            .lines()
            .zip(golden.lines())
            .position(|(a, b)| a != b)
            .unwrap_or_else(|| svg.lines().count().min(golden.lines().count()));
        panic!(
            "the twin paradox diagram no longer matches tests/golden/twin_paradox.svg; first difference at line {}:\n  now:    {}\n  golden: {}",
            line + 1,
            svg.lines().nth(line).unwrap_or("<end>"),
            golden.lines().nth(line).unwrap_or("<end>")
        );
    }
}

#[test]
fn the_twin_paradox_shows_the_right_ages() {
    let (distance, acceleration) = (4.2465, 1.0323);
    // one leg, rest to rest, from the core's closed form
    let leg = FlipBurn::new(distance, acceleration).unwrap();
    let diagram = twin_paradox_diagram(distance, acceleration).unwrap();
    let svg = diagram.render();
    let elements = parse(&svg).unwrap();
    let subtitle = &find(&elements, "subtitle", "text")[0].text;
    assert_eq!(
        subtitle,
        &format!(
            "At the reunion the home clock reads {:.2} yr and the ship clock {:.2} yr",
            2.0 * leg.t,
            2.0 * leg.tau
        )
    );
    assert!(subtitle.contains("11.75 yr") && subtitle.contains("7.08 yr"));

    let layout = diagram.layout();
    // the two named events are where the closed form puts them
    let events: Vec<(f64, f64)> = find(&elements, "event", "circle")
        .iter()
        .map(|c| c.centre())
        .collect();
    let near =
        |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).abs() < 0.011 && (a.1 - b.1).abs() < 0.011;
    assert_eq!(events.len(), 2);
    assert!(
        near(events[0], layout.to_pixels(distance, leg.t)),
        "turnaround"
    );
    assert!(
        near(events[1], layout.to_pixels(0.0, 2.0 * leg.t)),
        "reunion"
    );

    // a tick a year on each clock: 0 to 11 at home, 0 to 7 aboard
    let ticks = find(&elements, "ticks", "circle");
    assert_eq!(ticks.len(), 12 + 8);
    // every flash is light: 45 degrees, eleven from home and seven from the ship
    let flashes = find(&elements, "signal", "line");
    assert_eq!(flashes.len(), 11 + 7);
    for flash in &flashes {
        let ((x1, y1), (x2, y2)) = flash.ends();
        assert!(((x2 - x1).abs() - (y2 - y1).abs()).abs() <= 0.022);
        assert!(y2 < y1, "light goes up the page");
    }
    // flashes from home lean right, flashes from the ship lean left and end on the home worldline
    let home_x = layout.to_pixels(0.0, 0.0).0;
    assert_eq!(
        flashes
            .iter()
            .filter(|f| f.ends().1 .0 > f.ends().0 .0)
            .count(),
        11
    );
    assert!(flashes
        .iter()
        .filter(|f| f.ends().1 .0 < f.ends().0 .0)
        .all(|f| (f.ends().1 .0 - home_x).abs() <= 0.011));
    // each flash starts on a tick of the clock that sent it
    let tick_at = |p: (f64, f64)| ticks.iter().any(|t| near(t.centre(), p));
    assert!(flashes.iter().all(|f| tick_at(f.ends().0)));
}

#[test]
fn the_rindler_chase_sorts_the_flashes_by_the_horizon() {
    let acceleration = 0.5;
    let diagram = rindler_chase_diagram(acceleration).unwrap();
    let elements = parse(&diagram.render()).unwrap();
    let layout = diagram.layout();
    let groups: Vec<&Element> = elements
        .iter()
        .filter(|e| e.name == "g" && e.attr("class") == Some("signal"))
        .collect();
    // three flashes sent before t = 1/a arrive, two sent after never do: two colours of group
    assert_eq!(groups.len(), 5);
    let caught = groups[0].attr("stroke").unwrap();
    assert_eq!(
        groups
            .iter()
            .filter(|g| g.attr("stroke") == Some(caught))
            .count(),
        3
    );
    // the horizon passes through (t = 1/a, x = 0) at 45 degrees
    let horizon = find(&elements, "horizon", "line")[0].ends();
    let last_chance = layout.to_pixels(0.0, 1.0 / acceleration);
    assert!(common::distance_to_segment(last_chance, horizon.0, horizon.1) < 0.02);
    // every flash that arrives ends on a dot on the ship, and stays below the horizon line
    let arrivals: Vec<(f64, f64)> = find(&elements, "event", "circle")
        .iter()
        .map(|c| c.centre())
        .collect();
    assert_eq!(arrivals.len(), 3 + 1, "three arrivals and the last chance");
    let lines = find(&elements, "signal", "line");
    assert_eq!(lines.len(), 5);
    for line in &lines[..3] {
        let end = line.ends().1;
        assert!(arrivals
            .iter()
            .any(|a| (a.0 - end.0).abs() < 0.011 && (a.1 - end.1).abs() < 0.011));
    }
    // the others run off the top of the plot
    for line in &lines[3..] {
        assert!((line.ends().1 .1 - layout.plot_top).abs() < 0.011);
    }
    assert!(diagram.render().contains("at or after t = 2.00 yr"));
}

#[test]
fn the_scenes_refuse_bad_numbers_and_survive_extreme_ones() {
    for (distance, acceleration, error) in [
        (0.0, 1.0, Error::OutOfRange),
        (-1.0, 1.0, Error::OutOfRange),
        (1.0, 0.0, Error::OutOfRange),
        (1.0, -2.0, Error::OutOfRange),
        // too slow to tell from Newton's trip
        (1e-20, 1e-20, Error::OutOfRange),
        (f64::NAN, 1.0, Error::NonFinite),
        (1.0, f64::INFINITY, Error::NonFinite),
    ] {
        assert_eq!(
            twin_paradox(distance, acceleration).err(),
            Some(error),
            "{distance} {acceleration}"
        );
    }
    for (acceleration, error) in [
        (0.0, Error::OutOfRange),
        (-1.0, Error::OutOfRange),
        (1e-200, Error::OutOfRange),
        (1e200, Error::OutOfRange),
        (f64::NAN, Error::NonFinite),
    ] {
        assert_eq!(
            rindler_chase(acceleration).err(),
            Some(error),
            "{acceleration}"
        );
    }
    // Any other pair of numbers gives a well-formed picture or a clean refusal.
    let mut drawn = 0;
    for distance in [1e-9, 1e-4, 0.3, 7.0, 1e3, 1e7, 1e12, 1e150, 1e300] {
        for acceleration in [1e-9, 1e-3, 0.1, 1.0, 50.0, 1e6, 1e12, 1e300] {
            if let Ok(svg) = twin_paradox(distance, acceleration) {
                parse(&svg).unwrap_or_else(|e| panic!("{distance} {acceleration}: {e}"));
                assert!(
                    !svg.contains("NaN") && !svg.contains("inf"),
                    "{distance} {acceleration}"
                );
                drawn += 1;
            }
        }
    }
    assert!(drawn >= 40, "only {drawn} drawn");
    for acceleration in [1e-90, 1e-6, 0.01, 3.0, 1e4, 1e90] {
        let svg = rindler_chase(acceleration).unwrap();
        parse(&svg).unwrap();
        assert!(!svg.contains("NaN") && !svg.contains("inf"));
    }
}

#[test]
fn themes_change_the_colours_and_nothing_else() {
    let strip = |svg: &str| -> String {
        // drop every #rrggbb
        let mut out = String::new();
        let mut rest = svg;
        while let Some(hash) = rest.find('#') {
            out.push_str(&rest[..hash]);
            rest = &rest[hash + 7..];
        }
        out.push_str(rest);
        out
    };
    let (light, dark) = (
        everything(Theme::Light).render(),
        everything(Theme::Dark).render(),
    );
    assert_ne!(light, dark);
    assert_eq!(strip(&light), strip(&dark));
    assert!(light.contains("fill=\"#fcfcfb\"") && dark.contains("fill=\"#1a1a19\""));
    // a colour given outright is the same in both
    assert!(light.contains("#c81e1e") && dark.contains("#c81e1e"));
}

#[test]
fn the_legend_names_each_line_once_and_sits_beside_the_plot() {
    let diagram = twin_paradox_diagram(4.2465, 1.0323).unwrap();
    let elements = parse(&diagram.render()).unwrap();
    let layout = diagram.layout();
    let names: Vec<&str> = find(&elements, "legend", "text")
        .iter()
        .map(|t| t.text.as_str())
        .collect();
    assert_eq!(
        names,
        [
            "Twin at home",
            "Travelling twin",
            "Home clock, every 1 yr",
            "Ship clock, every 1 yr",
            "Light from home, one flash per tick",
            "Light from the ship, one flash per tick",
        ]
    );
    // everything in the legend is to the right of the plot, and inside the image
    let plot_right = layout.plot_left + layout.plot_width;
    for element in elements
        .iter()
        .filter(|e| e.classes.iter().any(|c| c == "legend") && e.name != "g")
    {
        let x = match element.name.as_str() {
            "line" => element.num("x1"),
            "circle" => element.num("cx") - element.num("r"),
            _ => element.num("x"),
        };
        assert!(x > plot_right + 10.0 && x < layout.width);
    }
    // swatches carry the dash pattern of the line they stand for: two series, two patterns
    let swatches = find(&elements, "legend", "line");
    assert_eq!(swatches.len(), 4);
    assert_eq!(swatches[0].attr("stroke-dasharray"), None);
    assert_eq!(swatches[1].attr("stroke-dasharray"), Some("12 4"));
    assert_ne!(
        swatches[2].attr("stroke-dasharray"),
        swatches[3].attr("stroke-dasharray")
    );

    // no legend when it is switched off or nothing is labelled
    let mut quiet = twin_paradox_diagram(4.2465, 1.0323).unwrap();
    let narrow = parse(&quiet.legend(false).render()).unwrap();
    assert!(find(&narrow, "legend", "text").is_empty());
    assert!(narrow[0].num("width") < layout.width);
}

#[test]
fn axis_ticks_are_at_nice_steps_and_do_not_crowd() {
    let labels = |view: View, axis: &str| -> Vec<(f64, String)> {
        let diagram = Diagram::new(view).unwrap();
        let elements = parse(&diagram.render()).unwrap();
        find(&elements, axis, "text")
            .iter()
            .filter(|t| t.attr("class") != Some("axis-title"))
            .map(|t| {
                (
                    if axis == "axis x" {
                        t.num("x")
                    } else {
                        t.num("y")
                    },
                    t.text.clone(),
                )
            })
            .collect()
    };
    let texts =
        |ticks: &[(f64, String)]| -> Vec<String> { ticks.iter().map(|t| t.1.clone()).collect() };

    let view = View::equal_scales((-0.85, 5.1), (-0.47, 12.2), 50.0);
    assert_eq!(texts(&labels(view, "axis x")), ["0", "2", "4"]);
    assert_eq!(
        texts(&labels(view, "axis t")),
        ["0", "2", "4", "6", "8", "10", "12"]
    );

    let cases = [
        View::new((0.0, 1.0), (0.0, 1.0), 400.0, 300.0),
        View::new((-0.03, 0.07), (1e6, 3.5e6), 640.0, 480.0),
        View::new((-1234.5, 987.6), (-2e-5, 9e-5), 900.0, 200.0),
        View::new((99_990.0, 100_010.0), (0.0, 1e12), 300.0, 700.0),
        View::new((0.0, 3.0), (0.0, 7.0), 60.0, 40.0),
    ];
    for view in cases {
        for axis in ["axis x", "axis t"] {
            let ticks = labels(view, axis);
            assert!(!ticks.is_empty(), "{view:?} {axis}");
            // read as numbers, the labels step by 1, 2 or 5 times a power of ten
            let values: Vec<f64> = ticks
                .iter()
                .map(|t| t.1.replace('\u{2212}', "-").parse().unwrap())
                .collect();
            if values.len() > 1 {
                let step = values[1] - values[0];
                assert!(step > 0.0);
                for pair in values.windows(2) {
                    assert!((pair[1] - pair[0] - step).abs() < 1e-9 * step.max(1.0));
                }
                let mantissa = step / 10f64.powf(step.log10().floor());
                assert!(
                    [1.0, 2.0, 5.0, 10.0]
                        .iter()
                        .any(|m| (mantissa - m).abs() < 1e-6),
                    "{view:?} {axis}: step {step}"
                );
                // and neighbours are far enough apart on the page to be read
                let gaps: Vec<f64> = ticks.windows(2).map(|p| (p[1].0 - p[0].0).abs()).collect();
                let least = if axis == "axis x" { 30.0 } else { 14.0 };
                assert!(
                    gaps.iter().all(|&g| g >= least),
                    "{view:?} {axis}: {gaps:?}"
                );
            }
            // no label is long: a power of ten is taken out of large and small numbers
            assert!(
                ticks.iter().all(|t| t.1.chars().count() <= 9),
                "{view:?} {axis}: {:?}",
                texts(&ticks)
            );
        }
    }
    // the power of ten taken out is shown beside the axis title
    let mut diagram = Diagram::new(View::new((0.0, 1.0), (1e6, 3.5e6), 300.0, 300.0)).unwrap();
    diagram.axis_titles("x", "t (s)");
    let elements = parse(&diagram.render()).unwrap();
    let title = find(&elements, "axis t", "text")
        .into_iter()
        .find(|t| t.attr("class") == Some("axis-title"))
        .unwrap();
    assert_eq!(title.text, "t (s) ×10");
    assert!(elements.iter().any(|e| e.name == "tspan" && e.text == "6"));
}
