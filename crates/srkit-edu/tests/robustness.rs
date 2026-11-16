//! Robustness: whatever numbers go in, every call returns a diagram or an error, and never panics.
mod common;

use common::parse;
use srkit_core::{
    kinematics::Event,
    segment::{Null, Segment},
    worldline::WorldlineBuf,
    Vec3,
};
use srkit_edu::{
    rindler_chase, twin_paradox, Colour, Cone, Dash, Diagram, Frame, Placement, Style, Theme,
    TickNumbers, View,
};

/// Numbers chosen to break things: not numbers, infinities, zeros of both signs, the smallest and largest
/// doubles, and a few ordinary ones.
const WILD: [f64; 17] = [
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    0.0,
    -0.0,
    5e-324,
    f64::MIN_POSITIVE,
    -1e-150,
    1e-9,
    -0.5,
    1.0,
    3.75,
    -1e6,
    1e15,
    1e150,
    -1e300,
    f64::MAX,
];

/// A repeatable pick from [`WILD`], so that the test needs no dependency and fails the same way twice.
struct Picks(u64);

impl Picks {
    fn index(&mut self, len: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) as usize) % len
    }
    fn wild(&mut self) -> f64 {
        WILD[self.index(WILD.len())]
    }
    /// Mostly ordinary numbers, sometimes a wild one: a diagram has to exist before much can go wrong in it.
    fn number(&mut self) -> f64 {
        if self.index(4) == 0 {
            self.wild()
        } else {
            [-7.0, -1.5, -0.25, 0.0, 0.3, 1.0, 2.5, 9.0][self.index(8)]
        }
    }
    fn vector(&mut self) -> Vec3 {
        Vec3::new(self.number(), self.number(), self.number())
    }
    fn event(&mut self) -> Event {
        Event::new(self.number(), self.vector())
    }
    fn style(&mut self) -> Style {
        let colour = match self.index(5) {
            0 => Colour::Series(self.index(256) as u8),
            1 => Colour::Light,
            2 => Colour::Ink,
            3 => Colour::Muted,
            _ => Colour::Rgb(self.index(256) as u8, 0, 255),
        };
        let dash = [
            Dash::Solid,
            Dash::LongDash,
            Dash::Dashed,
            Dash::DashDot,
            Dash::DashDotDot,
            Dash::Dotted,
        ][self.index(6)];
        let width = if self.index(6) == 0 { self.wild() } else { 1.5 };
        Style::new(colour, width, dash)
    }
    fn segment(&mut self) -> Segment {
        let (t0, t1, p) = (self.number(), self.number(), self.vector());
        match self.index(6) {
            0 => Segment::rest(t0, t1, p),
            1 => Segment::inertial(t0, t1, p, self.vector()),
            2 => Segment::hyperbolic(t0, t1, p, self.vector(), self.number(), self.number()),
            3 => Segment::circular(
                t0,
                t1,
                p,
                self.number(),
                self.number(),
                self.vector(),
                self.vector(),
                self.number(),
            ),
            4 => Segment::null(t0, t1, p, self.vector()),
            _ => Segment::Null(Null::between(t0, p, self.vector())),
        }
    }
}

#[test]
fn no_numbers_make_a_diagram_panic() {
    let mut picks = Picks(7);
    let (mut built, mut added, mut refused) = (0, 0, 0);
    for round in 0..1500 {
        // Ranges the right way round more often than not, so that some diagrams exist to be abused.
        let range = |picks: &mut Picks| {
            let (a, b) = (picks.number(), picks.number());
            if picks.index(5) == 0 {
                (a, b)
            } else {
                (a.min(b), a.max(b))
            }
        };
        let (x, t) = (range(&mut picks), range(&mut picks));
        let sizes = [1.0, 90.0, 640.0, picks.wild()];
        let mut view = View::new(x, t, sizes[picks.index(4)], sizes[picks.index(4)]);
        if picks.index(3) == 0 {
            view = view.axis(picks.vector(), picks.vector());
        }
        if picks.index(3) == 0 {
            view = view.frame(Frame::new(picks.event(), picks.number()));
        }
        let Ok(mut diagram) = Diagram::new(view) else {
            continue;
        };
        built += 1;
        diagram
            .theme(if round % 2 == 0 {
                Theme::Light
            } else {
                Theme::Dark
            })
            .title("wild")
            .axis_titles("x", "t")
            .legend(round % 3 != 0)
            .grid(round % 5 != 0);
        for _ in 0..12 {
            let style = picks.style();
            let (a, b) = (picks.event(), picks.event());
            let frame = Frame::new(picks.event(), picks.number());
            let (segment, next) = (picks.segment(), picks.segment());
            let worldline = WorldlineBuf::new([segment, next], picks.number(), picks.number());
            let outcome = match picks.index(12) {
                0 => diagram.event(a, "event", style).err(),
                1 => diagram
                    .event_placed(a, "event", Placement::Below, style)
                    .err(),
                2 => diagram.light_cone(a, Cone::Both, "cone", style).err(),
                3 => diagram.signal(a, b, "signal", style).err(),
                4 => diagram
                    .simultaneity(frame, picks.number(), "now", style)
                    .err(),
                5 => diagram
                    .constant_position(frame, picks.number(), "here", style)
                    .err(),
                6 => diagram.frame_axes(frame, "axes", style).err(),
                7 => diagram
                    .frame_grid(frame, picks.number(), "grid", style)
                    .err(),
                8 => diagram.rindler_horizon(&segment, "horizon", style).err(),
                9 => match &worldline {
                    Ok(w) => diagram.worldline(&w.view(), "worldline", style).err(),
                    Err(_) => continue,
                },
                10 => match &worldline {
                    Ok(w) => diagram
                        .proper_time_ticks(
                            &w.view(),
                            picks.number(),
                            TickNumbers::Right,
                            "ticks",
                            style,
                        )
                        .err(),
                    Err(_) => continue,
                },
                _ => match &worldline {
                    Ok(w) => diagram.seen(a, &w.view(), "seen", style).err(),
                    Err(_) => continue,
                },
            };
            match outcome {
                None => added += 1,
                Some(_) => refused += 1,
            }
        }
        // Whatever was accepted renders to a well-formed document with finite numbers in it.
        let svg = diagram.render();
        parse(&svg).unwrap_or_else(|e| panic!("round {round}: {e}"));
        assert!(
            !svg.contains("NaN") && !svg.contains("inf"),
            "round {round}"
        );
        let layout = diagram.layout();
        assert!(layout.width.is_finite() && layout.height.is_finite());
    }
    // the test did exercise both the accepting and the refusing paths
    assert!(built > 300, "only {built} diagrams built");
    assert!(
        added > 500 && refused > 500,
        "{added} added, {refused} refused"
    );
}

#[test]
fn no_numbers_make_a_scene_panic() {
    for &distance in &WILD {
        for &acceleration in &WILD {
            if let Ok(svg) = twin_paradox(distance, acceleration) {
                parse(&svg).unwrap_or_else(|e| panic!("{distance} {acceleration}: {e}"));
                assert!(distance > 0.0 && acceleration > 0.0);
            }
        }
        if let Ok(svg) = rindler_chase(distance) {
            parse(&svg).unwrap_or_else(|e| panic!("{distance}: {e}"));
            assert!(distance > 0.0);
        }
    }
}
