//! The diagram builder: what is drawn, and how it becomes SVG.

use srkit_core::{
    kinematics::Event, rindler::RindlerHorizon, segment::Segment, worldline::Worldline,
};

use crate::{
    axis::{Readings, Step, Ticks},
    geometry::{Pt, Rect},
    sample::{natural_range, sample, visible_range, Window},
    style::{hex, Colour, Style, Theme},
    svg::{fixed, text_width, Svg},
    view::{check_frame, Chart, Frame, View},
    Error,
};

/// The most ticks or grid lines one call may add.
pub const MAX_MARKS: usize = 2000;

/// The relative tolerance of [`Diagram::signal`]: the time between the two events may differ from their
/// distance by this fraction of the larger of the two. A further 10⁻¹² of the size of the coordinates
/// themselves is allowed for rounding, which matters only for events far from the origin and close together.
pub const NULL_TOLERANCE: f64 = 1e-9;

/// System fonts only, so nothing is fetched. Family names with a space need no quotes in CSS.
const FONT: &str = "system-ui, -apple-system, Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif";
/// Space kept clear round the edge of the image.
const PAD: f64 = 14.0;
const TITLE_SIZE: f64 = 15.0;
const SUBTITLE_SIZE: f64 = 12.0;
const AXIS_SIZE: f64 = 12.0;
const TICK_SIZE: f64 = 11.0;
const MARK_SIZE: f64 = 11.0;
const LEGEND_SIZE: f64 = 12.0;
const TICK_LENGTH: f64 = 5.0;
/// From the bottom of the plot down to the baseline of the tick labels of the position axis: clear of the
/// lowest label of the time axis, which sits beside the corner.
const X_LABEL_DROP: f64 = TICK_LENGTH + 5.0 + TICK_SIZE - 2.0;
/// Least distance between ticks on the position axis and on the time axis.
const X_TICK_GAP: f64 = 56.0;
const T_TICK_GAP: f64 = 34.0;
/// Least clear space between two tick labels on the position axis.
const X_LABEL_GAP: f64 = 14.0;
const LEGEND_GAP: f64 = 22.0;
const LEGEND_SWATCH: f64 = 34.0;
const LEGEND_ROW: f64 = 21.0;
const EVENT_RADIUS: f64 = 4.0;
const TICK_RADIUS: f64 = 3.0;

/// Which half of a light cone to draw.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cone {
    /// The events whose light reaches the apex.
    Past,
    /// The events the apex's light reaches.
    Future,
    /// Both halves.
    Both,
}

/// Where a label sits beside its dot.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Placement {
    /// To the right, or to the left where the right would run past the plot.
    #[default]
    Auto,
    /// To the right of the dot.
    Right,
    /// To the left of the dot.
    Left,
    /// Above the dot.
    Above,
    /// Below the dot.
    Below,
}

/// Whether proper-time ticks carry their clock reading, and on which side.
///
/// A reading is left out where it would run into a dot, into a label already placed, or off the plot, so
/// closely spaced ticks keep their dots and lose some of their numbers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TickNumbers {
    /// Dots only.
    #[default]
    None,
    /// The proper time to the left of each tick.
    Left,
    /// The proper time to the right of each tick.
    Right,
}

/// Where the plot sits in the image, for placing things over a rendered diagram.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Layout {
    /// Width of the whole image in pixels.
    pub width: f64,
    /// Height of the whole image in pixels.
    pub height: f64,
    /// Left edge of the plot area.
    pub plot_left: f64,
    /// Top edge of the plot area.
    pub plot_top: f64,
    /// Width of the plot area.
    pub plot_width: f64,
    /// Height of the plot area.
    pub plot_height: f64,
    /// The position range across the plot.
    pub x: (f64, f64),
    /// The time range up the plot.
    pub t: (f64, f64),
}

impl Layout {
    /// The pixel position of the point with the given coordinates in the drawing frame. Pixel `y` grows
    /// downwards, so later times have smaller `y`.
    pub fn to_pixels(&self, x: f64, t: f64) -> (f64, f64) {
        (
            self.plot_left + (x - self.x.0) * (self.plot_width / (self.x.1 - self.x.0)),
            self.plot_top + (self.t.1 - t) * (self.plot_height / (self.t.1 - self.t.0)),
        )
    }
    fn px(&self, p: Pt) -> (f64, f64) {
        self.to_pixels(p.x, p.t)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Swatch {
    None,
    Line,
    Dot,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Line {
    a: Pt,
    b: Pt,
    /// Draw an arrowhead half way along, pointing from `a` to `b`.
    arrow: bool,
}

#[derive(Clone, Debug, PartialEq)]
struct Dot {
    p: Pt,
    text: String,
    placement: Placement,
    radius: f64,
    /// A ring in the item's colour instead of a filled disc.
    hollow: bool,
    /// Whether the text is a clock reading (drawn in the quieter ink) or a name.
    quiet: bool,
}

/// One thing added to the diagram, already reduced to geometry in the drawing frame and clipped to the view.
#[derive(Clone, Debug, PartialEq)]
struct Item {
    /// Class of the group that holds the lines.
    class: &'static str,
    /// Class of the group that holds the dots.
    dot_class: &'static str,
    label: String,
    style: Style,
    swatch: Swatch,
    paths: Vec<Vec<Pt>>,
    lines: Vec<Line>,
    dots: Vec<Dot>,
}

impl Item {
    fn new(class: &'static str, label: &str, style: Style, swatch: Swatch) -> Self {
        Self {
            class,
            dot_class: class,
            label: label.to_string(),
            style,
            swatch,
            paths: Vec::new(),
            lines: Vec::new(),
            dots: Vec::new(),
        }
    }
}

/// A 1+1D Minkowski diagram: time up, one spatial axis across.
///
/// Build one from a [`View`], add worldlines, events and light, then call [`Diagram::render`] for the SVG.
/// Everything is given in the lab frame, the frame the worldlines are stored in; the diagram transforms it
/// into the frame of its view and projects it onto the view's axis. Each thing takes a label and a [`Style`].
/// A labelled line or set of ticks gets an entry in the legend (identical entries are merged, and an empty
/// label gives none); an event's label is written beside its dot.
///
/// The same calls in the same order always give the same bytes.
#[derive(Clone, Debug, PartialEq)]
pub struct Diagram {
    view: View,
    chart: Chart,
    rect: Rect,
    theme: Theme,
    title: String,
    subtitle: String,
    x_title: String,
    t_title: String,
    legend: bool,
    grid: bool,
    items: Vec<Item>,
}

fn check_style(style: Style) -> Result<(), Error> {
    if style.width.is_finite() && style.width > 0.0 && style.width <= 64.0 {
        Ok(())
    } else {
        Err(Error::InvalidStyle)
    }
}

fn check_event(event: Event) -> Result<(), Error> {
    if event.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFinite)
    }
}

impl Diagram {
    /// An empty diagram of `view`, in the light theme, with a legend and a background grid.
    pub fn new(view: View) -> Result<Self, Error> {
        let view = view.validated()?;
        Ok(Self {
            view,
            chart: Chart::new(view.origin, view.direction, view.frame),
            rect: Rect {
                x0: view.x.0,
                x1: view.x.1,
                t0: view.t.0,
                t1: view.t.1,
            },
            theme: Theme::Light,
            title: String::new(),
            subtitle: String::new(),
            x_title: String::new(),
            t_title: String::new(),
            legend: true,
            grid: true,
            items: Vec::new(),
        })
    }

    /// Choose the colour scheme.
    pub fn theme(&mut self, theme: Theme) -> &mut Self {
        self.theme = theme;
        self
    }
    /// Set the title, written above the plot. Empty for none.
    pub fn title(&mut self, title: &str) -> &mut Self {
        self.title = title.to_string();
        self
    }
    /// Set a second, smaller line under the title. Empty for none.
    pub fn subtitle(&mut self, subtitle: &str) -> &mut Self {
        self.subtitle = subtitle.to_string();
        self
    }
    /// Set the axis titles. The diagram knows no units, so name them here: `"x (ly)"`, `"t (yr)"`.
    pub fn axis_titles(&mut self, x: &str, t: &str) -> &mut Self {
        self.x_title = x.to_string();
        self.t_title = t.to_string();
        self
    }
    /// Show or hide the legend. It sits to the right of the plot, never over it.
    pub fn legend(&mut self, show: bool) -> &mut Self {
        self.legend = show;
        self
    }
    /// Show or hide the faint grid lines at the axis ticks.
    pub fn grid(&mut self, show: bool) -> &mut Self {
        self.grid = show;
        self
    }

    /// The coordinates `(x, t)` of a lab event in the frame the diagram is drawn in.
    pub fn coordinates(&self, event: Event) -> (f64, f64) {
        let p = self.chart.coords(event);
        (p.x, p.t)
    }

    fn window(&self) -> Window {
        Window {
            chart: self.chart,
            rect: self.rect,
            sx: self.view.width / (self.rect.x1 - self.rect.x0),
            st: self.view.height / (self.rect.t1 - self.rect.t0),
        }
    }

    fn point(&self, event: Event) -> Result<Pt, Error> {
        check_event(event)?;
        let p = self.chart.coords(event);
        if p.is_finite() {
            Ok(p)
        } else {
            Err(Error::NonFinite)
        }
    }

    /// Add a worldline over the whole of its own range: from the start of its first segment to the end of
    /// its last (see [`Diagram::worldline_between`] to draw less).
    ///
    /// Straight segments (rest, inertial, light-like) take two points. Curved ones are cut until the drawn
    /// line is within [`crate::CURVE_TOLERANCE`] pixels of the true curve. Nothing is evaluated outside the
    /// worldline's own range or outside the view's time range.
    pub fn worldline(
        &mut self,
        worldline: &Worldline<'_>,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        self.worldline_between(worldline, f64::NEG_INFINITY, f64::INFINITY, label, style)
    }

    /// Add the part of a worldline between lab times `from` and `to`. The part outside the worldline's own
    /// range is left out.
    pub fn worldline_between(
        &mut self,
        worldline: &Worldline<'_>,
        from: f64,
        to: f64,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        if from.is_nan() || to.is_nan() {
            return Err(Error::NonFinite);
        }
        let window = self.window();
        let mut item = Item::new("worldline", label, style, Swatch::Line);
        if let Some((a, b)) = visible_range(worldline, &window, from, to)? {
            item.paths = self.rect.clip_polyline(&sample(worldline, &window, a, b)?);
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add an event: a dot with its label beside it. An event outside the view is not drawn.
    pub fn event(&mut self, event: Event, label: &str, style: Style) -> Result<&mut Self, Error> {
        self.event_placed(event, label, Placement::Auto, style)
    }

    /// Add an event and say where its label goes.
    pub fn event_placed(
        &mut self,
        event: Event,
        label: &str,
        placement: Placement,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        let p = self.point(event)?;
        let mut item = Item::new("event", "", style, Swatch::None);
        if self.rect.contains(p) {
            item.dots.push(Dot {
                p,
                text: label.to_string(),
                placement,
                radius: EVENT_RADIUS,
                hollow: false,
                quiet: false,
            });
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add the light cone of an event: the two light rays into its past, out to its future, or both, as far
    /// as the edge of the view. With equal scales they are at 45° in every frame.
    pub fn light_cone(
        &mut self,
        apex: Event,
        cone: Cone,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        let p = self.point(apex)?;
        let mut item = Item::new("cone", label, style, Swatch::Line);
        let halves: &[f64] = match cone {
            Cone::Past => &[-1.0],
            Cone::Future => &[1.0],
            Cone::Both => &[-1.0, 1.0],
        };
        for &t in halves {
            for x in [-1.0, 1.0] {
                if let Some((a, b)) = self.rect.clip_ray(p, Pt { x, t }) {
                    item.lines.push(Line { a, b, arrow: false });
                }
            }
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add a light signal that leaves the event `from` and arrives at the event `to`, with an arrowhead half
    /// way along.
    ///
    /// The time from `from` to `to` must equal the distance between them to within [`NULL_TOLERANCE`];
    /// otherwise the two are not joined by light and the result is [`Error::NotLightLike`]. That includes
    /// a pair given in the wrong order, since light does not arrive before it leaves.
    pub fn signal(
        &mut self,
        from: Event,
        to: Event,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        let (a, b) = (self.point(from)?, self.point(to)?);
        let (dt, distance) = (to.t - from.t, to.x.dist(from.x));
        let size = from.t.abs() + to.t.abs() + from.x.norm() + to.x.norm();
        let tolerance = NULL_TOLERANCE * dt.abs().max(distance) + 1e-12 * size;
        if !((dt - distance).abs() <= tolerance) {
            return Err(Error::NotLightLike);
        }
        let mut item = Item::new("signal", label, style, Swatch::Line);
        if let Some((a, b)) = self.rect.clip_segment(a, b) {
            item.lines.push(Line { a, b, arrow: true });
        }
        self.items.push(item);
        Ok(self)
    }

    /// The chart of another frame along this diagram's axis.
    fn frame_chart(&self, frame: Frame) -> Result<Chart, Error> {
        check_frame(frame)?;
        Ok(Chart::new(
            self.view.origin,
            self.view.direction,
            Some(frame),
        ))
    }

    /// The line on which one coordinate of `frame` is constant: its time (`time = true`) or its position.
    fn frame_line(&self, frame: &Chart, time: bool, value: f64) -> Result<Option<Line>, Error> {
        // A point of the line in the lab, then in the drawing frame.
        let on_line = frame.event(if time {
            Pt {
                x: frame.pivot_x,
                t: value,
            }
        } else {
            Pt {
                x: value,
                t: frame.pivot.t,
            }
        });
        let p = self.chart.coords(on_line);
        // Seen from the drawing frame the other frame moves at the difference of the rapidities: its lines of
        // simultaneity rise by that speed per unit of position, its lines of position lean by it per unit of time.
        let speed = libm::tanh(frame.rapidity - self.chart.rapidity);
        if !p.is_finite() || !speed.is_finite() {
            return Err(Error::NonFinite);
        }
        let d = if time {
            Pt { x: 1.0, t: speed }
        } else {
            Pt { x: speed, t: 1.0 }
        };
        Ok(self
            .rect
            .clip_line(p, d)
            .map(|(a, b)| Line { a, b, arrow: false }))
    }

    /// Add a line of simultaneity of `frame`: the events that share the time `t` in that frame.
    ///
    /// Drawn in another frame it is tilted: it rises by the frame's relative speed for each unit of position.
    pub fn simultaneity(
        &mut self,
        frame: Frame,
        t: f64,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        self.frame_lines(frame, &[(true, t)], "simultaneity", label, style)
    }

    /// Add a line of constant position of `frame`: the worldline of a point at rest at position `x` in that
    /// frame.
    pub fn constant_position(
        &mut self,
        frame: Frame,
        x: f64,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        self.frame_lines(frame, &[(false, x)], "position", label, style)
    }

    /// Add the two axes of `frame` through its origin event: its time axis and its line of simultaneity
    /// there. These are the tilted axes of a moving observer.
    pub fn frame_axes(
        &mut self,
        frame: Frame,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        let chart = self.frame_chart(frame)?;
        self.frame_lines(
            frame,
            &[(true, chart.pivot.t), (false, chart.pivot_x)],
            "frame-axes",
            label,
            style,
        )
    }

    fn frame_lines(
        &mut self,
        frame: Frame,
        lines: &[(bool, f64)],
        class: &'static str,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        let chart = self.frame_chart(frame)?;
        let mut item = Item::new(class, label, style, Swatch::Line);
        for &(time, value) in lines {
            if !value.is_finite() {
                return Err(Error::NonFinite);
            }
            item.lines.extend(self.frame_line(&chart, time, value)?);
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add the coordinate grid of `frame`: its lines of simultaneity and of constant position at every
    /// multiple of `step` that crosses the view. [`Error::TooManyMarks`] beyond [`MAX_MARKS`] lines.
    pub fn frame_grid(
        &mut self,
        frame: Frame,
        step: f64,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        if !step.is_finite() {
            return Err(Error::NonFinite);
        }
        if !(step > 0.0) {
            return Err(Error::OutOfRange);
        }
        let chart = self.frame_chart(frame)?;
        // The frame's coordinates of the corners of the view bound those of everything in it.
        let r = self.rect;
        let (mut t_range, mut x_range) = (
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::INFINITY, f64::NEG_INFINITY),
        );
        for (x, t) in [(r.x0, r.t0), (r.x0, r.t1), (r.x1, r.t0), (r.x1, r.t1)] {
            let c = chart.coords(self.chart.event(Pt { x, t }));
            if !c.is_finite() {
                return Err(Error::NonFinite);
            }
            t_range = (t_range.0.min(c.t), t_range.1.max(c.t));
            x_range = (x_range.0.min(c.x), x_range.1.max(c.x));
        }
        let multiples = |(lo, hi): (f64, f64)| {
            let first = libm::ceil(lo / step);
            let count = libm::floor(hi / step) - first + 1.0;
            (first, count.max(0.0))
        };
        let (first_t, count_t) = multiples(t_range);
        let (first_x, count_x) = multiples(x_range);
        if !(count_t + count_x <= MAX_MARKS as f64) {
            return Err(Error::TooManyMarks);
        }
        let mut item = Item::new("frame-grid", label, style, Swatch::Line);
        for (time, first, count) in [(true, first_t, count_t), (false, first_x, count_x)] {
            for i in 0..count as usize {
                let value = (first + i as f64) * step;
                item.lines.extend(self.frame_line(&chart, time, value)?);
            }
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add ticks along a worldline at every multiple of `dtau` of its proper time, as dots, with the clock
    /// reading beside each if `numbers` asks for it.
    ///
    /// The ticks are placed with `Worldline::t_of_tau`, so on a clock moving at Lorentz factor γ they are
    /// `γ dtau` apart in coordinate time. Only the part of the worldline that [`Diagram::worldline`] draws is
    /// ticked. [`Error::TooManyMarks`] beyond [`MAX_MARKS`] ticks.
    pub fn proper_time_ticks(
        &mut self,
        worldline: &Worldline<'_>,
        dtau: f64,
        numbers: TickNumbers,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        if !dtau.is_finite() {
            return Err(Error::NonFinite);
        }
        if !(dtau > 0.0) {
            return Err(Error::OutOfRange);
        }
        let window = self.window();
        let mut item = Item::new("ticks", label, style, Swatch::Dot);
        if let Some((a, b)) = visible_range(worldline, &window, f64::NEG_INFINITY, f64::INFINITY)? {
            let (tau_a, tau_b) = (worldline.tau(a), worldline.tau(b));
            if !tau_a.is_finite() || !tau_b.is_finite() {
                return Err(Error::NonFinite);
            }
            let first = libm::ceil(tau_a / dtau - 1e-9);
            let count = libm::floor(tau_b / dtau + 1e-9) - first + 1.0;
            if !(count <= MAX_MARKS as f64) {
                return Err(Error::TooManyMarks);
            }
            // `as` saturates, and the count is at most MAX_MARKS.
            let readings = Readings::new(dtau, first as i64, (first + count - 1.0) as i64);
            let placement = match numbers {
                TickNumbers::None => None,
                TickNumbers::Left => Some(Placement::Left),
                TickNumbers::Right => Some(Placement::Right),
            };
            for i in 0..count.max(0.0) as usize {
                let tau = (first + i as f64) * dtau;
                // Rounding may put the first or last tick a hair outside the stretch that is drawn.
                let t = worldline.t_of_tau(tau).max(a).min(b);
                let p = self.point(Event::new(t, worldline.at(t).x))?;
                if self.rect.contains(p) {
                    item.dots.push(Dot {
                        p,
                        text: placement
                            .map_or(String::new(), |_| readings.text((first + i as f64) as i64)),
                        placement: placement.unwrap_or_default(),
                        radius: TICK_RADIUS,
                        hollow: false,
                        quiet: true,
                    });
                }
            }
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add the Rindler horizon of a segment of constant proper acceleration, continued for ever: the line
    /// beyond which no signal reaches the ship while it keeps accelerating.
    ///
    /// For a ship accelerating along the diagram's axis it is the light ray through the pivot of the
    /// hyperbola, running the way the ship is pushed. A segment of any other kind has no horizon and gives
    /// [`Error::NoHorizon`].
    pub fn rindler_horizon(
        &mut self,
        segment: &Segment,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        // The horizon is the null plane (x - pivot.x)·n = t - pivot.t for a unit vector n.
        let (pivot, n) = match segment {
            Segment::Hyperbolic(g) => {
                let horizon = RindlerHorizon::of(g);
                (horizon.pivot, horizon.n)
            }
            Segment::Steered(g) => (
                Event::new(g.t0 - g.nt / g.alpha, g.p0 - g.ns / g.alpha),
                (g.u0 + g.ns) / g.plus,
            ),
            _ => return Err(Error::NoHorizon),
        };
        if !pivot.is_finite() || !n.is_finite() {
            return Err(Error::NonFinite);
        }
        // Two events on the axis that lie on the horizon.
        let (origin, dir) = (self.view.origin, self.view.direction);
        let on_axis = |s: f64| {
            let x = origin + dir * s;
            Event::new(pivot.t + (x - pivot.x).dot(n), x)
        };
        let s = (pivot.x - origin).dot(dir);
        let span = self.rect.x1 - self.rect.x0;
        let (a, b) = (self.point(on_axis(s))?, self.point(on_axis(s + span))?);
        let mut item = Item::new("horizon", label, style, Swatch::Line);
        let d = Pt {
            x: b.x - a.x,
            t: b.t - a.t,
        };
        if let Some((a, b)) = self.rect.clip_line(a, d) {
            item.lines.push(Line { a, b, arrow: false });
        }
        self.items.push(item);
        Ok(self)
    }

    /// Add the "what is seen" construction: the light ray that reaches `observer` from `source`, from the
    /// event where it left the source's worldline (found with `Worldline::send`, and marked with a ring) to
    /// the observer, with an arrowhead half way along.
    ///
    /// [`Error::NothingSeen`] if no light from the drawn part of the source reaches the observer.
    pub fn seen(
        &mut self,
        observer: Event,
        source: &Worldline<'_>,
        label: &str,
        style: Style,
    ) -> Result<&mut Self, Error> {
        check_style(style)?;
        let b = self.point(observer)?;
        let emitted = source.send(observer);
        let (lo, hi) = natural_range(source);
        if !(emitted.is_finite() && emitted >= lo && emitted <= hi) {
            return Err(Error::NothingSeen);
        }
        let a = self.point(Event::new(emitted, source.at(emitted).x))?;
        let mut item = Item::new("seen", label, style, Swatch::Line);
        item.dot_class = "seen-source";
        if let Some((a, b)) = self.rect.clip_segment(a, b) {
            item.lines.push(Line { a, b, arrow: true });
        }
        if self.rect.contains(a) {
            item.dots.push(Dot {
                p: a,
                text: String::new(),
                placement: Placement::Auto,
                radius: EVENT_RADIUS,
                hollow: true,
                quiet: true,
            });
        }
        self.items.push(item);
        Ok(self)
    }

    /// Where the plot sits in the rendered image.
    pub fn layout(&self) -> Layout {
        self.furniture().layout
    }

    /// The diagram as one self-contained SVG document: no scripts, no external fonts or references, every
    /// label escaped and every number written with two decimals.
    pub fn render(&self) -> String {
        let furniture = self.furniture();
        let mut svg = Svg::new();
        self.write(&mut svg, &furniture);
        svg.finish()
    }
}

#[derive(Clone, Debug, PartialEq)]
struct LegendEntry {
    label: String,
    style: Style,
    swatch: Swatch,
}

/// Everything round the plot, worked out before any of it is written.
struct Furniture {
    layout: Layout,
    x_ticks: Ticks,
    t_ticks: Ticks,
    legend: Vec<LegendEntry>,
    title_base: f64,
    subtitle_base: f64,
}

fn tick_labels(ticks: &Ticks) -> Vec<String> {
    ticks.indices().map(|k| ticks.label(k)).collect()
}

fn widest(labels: &[String], size: f64) -> f64 {
    labels
        .iter()
        .map(|l| text_width(l, size))
        .fold(0.0, f64::max)
}

impl Diagram {
    fn legend_entries(&self) -> Vec<LegendEntry> {
        let mut entries: Vec<LegendEntry> = Vec::new();
        if !self.legend {
            return entries;
        }
        for item in &self.items {
            if item.label.is_empty() || item.swatch == Swatch::None {
                continue;
            }
            let entry = LegendEntry {
                label: item.label.clone(),
                style: item.style,
                swatch: item.swatch,
            };
            if !entries.contains(&entry) {
                entries.push(entry);
            }
        }
        entries
    }

    fn furniture(&self) -> Furniture {
        let v = &self.view;
        let (x_span, t_span) = (v.x.1 - v.x.0, v.t.1 - v.t.0);
        let (sx, st) = (v.width / x_span, v.height / t_span);
        let equal = (sx - st).abs() <= 1e-9 * sx.max(st);

        // Tick steps. With equal scales both axes take the larger step, so the grid is made of squares.
        let mut x_step = Step::at_least(x_span * X_TICK_GAP / v.width);
        let mut t_step = Step::at_least(t_span * T_TICK_GAP / v.height);
        if equal {
            let larger = if x_step.value() >= t_step.value() {
                x_step
            } else {
                t_step
            };
            (x_step, t_step) = (larger, larger);
        }
        // Widen the position step until its labels have room side by side.
        let mut x_ticks = Ticks::new(v.x.0, v.x.1, x_step);
        for _ in 0..64 {
            let room = x_step.value() * sx - widest(&tick_labels(&x_ticks), TICK_SIZE);
            if x_ticks.count() <= 2 || room >= X_LABEL_GAP {
                break;
            }
            x_step = x_step.next();
            x_ticks = Ticks::new(v.x.0, v.x.1, x_step);
        }
        if equal && x_step.value() > t_step.value() {
            t_step = x_step;
        }
        let t_ticks = Ticks::new(v.t.0, v.t.1, t_step);
        let (x_labels, t_labels) = (tick_labels(&x_ticks), tick_labels(&t_ticks));

        // Left: the time axis title, its tick labels and the tick marks.
        let t_titled = !self.t_title.is_empty() || t_ticks.exponent != 0;
        let mut plot_left = PAD + widest(&t_labels, TICK_SIZE) + TICK_LENGTH + 4.0;
        if t_titled {
            plot_left += AXIS_SIZE + 8.0;
        }
        let first_x = x_labels.first().map_or(0.0, |l| text_width(l, TICK_SIZE));
        // The plot starts on a whole pixel.
        plot_left = plot_left.max(PAD + 0.5 * first_x).ceil();

        // Top: title and subtitle.
        let mut y = PAD;
        let mut title_base = 0.0;
        let mut subtitle_base = 0.0;
        if !self.title.is_empty() {
            title_base = y + TITLE_SIZE - 2.0;
            y = title_base + 8.0;
        }
        if !self.subtitle.is_empty() {
            subtitle_base = y + SUBTITLE_SIZE - 2.0;
            y = subtitle_base + 8.0;
        }
        if y > PAD {
            y += 6.0;
        }
        let plot_top = y.ceil();
        let (plot_right, plot_bottom) = (plot_left + v.width, plot_top + v.height);

        // Right: the legend, or room for the last tick label.
        let legend = self.legend_entries();
        let last_x = x_labels.last().map_or(0.0, |l| text_width(l, TICK_SIZE));
        let mut width = plot_right + PAD.max(0.5 * last_x + 4.0);
        if !legend.is_empty() {
            let text = legend
                .iter()
                .map(|e| text_width(&e.label, LEGEND_SIZE))
                .fold(0.0, f64::max);
            width = width.max(plot_right + LEGEND_GAP + LEGEND_SWATCH + 8.0 + text + PAD);
        }
        width = width
            .max(PAD + text_width(&self.title, TITLE_SIZE) * 1.06 + PAD)
            .max(PAD + text_width(&self.subtitle, SUBTITLE_SIZE) + PAD);

        // Bottom: tick labels and the position axis title.
        let x_titled = !self.x_title.is_empty() || x_ticks.exponent != 0;
        let mut height = plot_bottom + X_LABEL_DROP + PAD;
        if x_titled {
            height += AXIS_SIZE + 7.0;
        }
        let legend_bottom = plot_top + 4.0 + legend.len() as f64 * LEGEND_ROW;
        height = height.max(legend_bottom + PAD);

        Furniture {
            layout: Layout {
                width: width.ceil(),
                height: height.ceil(),
                plot_left,
                plot_top,
                plot_width: v.width,
                plot_height: v.height,
                x: v.x,
                t: v.t,
            },
            x_ticks,
            t_ticks,
            legend,
            title_base,
            subtitle_base,
        }
    }

    fn colour(&self, colour: Colour) -> String {
        hex(self.theme.rgb(colour))
    }

    /// The opening tag of a group that strokes its children in `style`.
    fn stroke_group(&self, svg: &mut Svg, class: &str, style: Style) {
        svg.open("g")
            .attr("class", class)
            .attr("fill", "none")
            .attr("stroke", &self.colour(style.colour))
            .num("stroke-width", style.width)
            .attr("stroke-linejoin", "round");
        if let Some(pattern) = style.dash.pattern() {
            svg.attr("stroke-dasharray", pattern);
        }
        svg.children();
    }

    fn tooltip(svg: &mut Svg, label: &str) {
        if !label.is_empty() {
            svg.open("title").content().text(label).close("title");
        }
    }

    fn write(&self, svg: &mut Svg, f: &Furniture) {
        let palette = self.theme.palette();
        let lay = &f.layout;
        let (surface, ink, secondary, muted, grid) = (
            hex(palette.surface),
            hex(palette.ink),
            hex(palette.secondary),
            hex(palette.muted),
            hex(palette.grid),
        );
        let (left, top) = (lay.plot_left, lay.plot_top);
        let (right, bottom) = (left + lay.plot_width, top + lay.plot_height);
        let name = if self.title.is_empty() {
            "Spacetime diagram"
        } else {
            &self.title
        };

        svg.open("svg")
            .attr("xmlns", "http://www.w3.org/2000/svg")
            .num("width", lay.width)
            .num("height", lay.height)
            .attr(
                "viewBox",
                &format!("0 0 {} {}", fixed(lay.width, 2), fixed(lay.height, 2)),
            )
            .attr("role", "img")
            .attr("aria-label", name)
            .attr("font-family", FONT)
            .children();
        svg.open("title").content().text(name).close("title");
        svg.open("desc")
            .content()
            .text(&self.description(f))
            .close("desc");
        svg.open("rect")
            .attr("class", "surface")
            .num("width", lay.width)
            .num("height", lay.height)
            .attr("fill", &surface)
            .empty();

        // The grid and the frame of the plot.
        let x_at = |k: i64| {
            lay.to_pixels(f.x_ticks.value(k), lay.t.0)
                .0
                .max(left)
                .min(right)
        };
        let y_at = |k: i64| {
            lay.to_pixels(lay.x.0, f.t_ticks.value(k))
                .1
                .max(top)
                .min(bottom)
        };
        if self.grid {
            svg.open("g")
                .attr("class", "grid")
                .attr("stroke", &grid)
                .num("stroke-width", 1.0)
                .children();
            for k in f.x_ticks.indices() {
                let x = x_at(k);
                svg.open("line")
                    .num("x1", x)
                    .num("y1", top)
                    .num("x2", x)
                    .num("y2", bottom)
                    .empty();
            }
            for k in f.t_ticks.indices() {
                let y = y_at(k);
                svg.open("line")
                    .num("x1", left)
                    .num("y1", y)
                    .num("x2", right)
                    .num("y2", y)
                    .empty();
            }
            svg.close("g");
        }
        svg.open("rect")
            .attr("class", "frame")
            .num("x", left)
            .num("y", top)
            .num("width", lay.plot_width)
            .num("height", lay.plot_height)
            .attr("fill", "none")
            .attr("stroke", &muted)
            .num("stroke-width", 1.0)
            .empty();

        // Lines first, then dots, then the text beside the dots, so that nothing is drawn over a label.
        for item in &self.items {
            if item.paths.is_empty() && item.lines.is_empty() {
                continue;
            }
            self.stroke_group(svg, item.class, item.style);
            Self::tooltip(svg, &item.label);
            for path in &item.paths {
                svg.open("polyline")
                    .points(path.iter().map(|&p| lay.px(p)))
                    .empty();
            }
            for line in &item.lines {
                let ((x1, y1), (x2, y2)) = (lay.px(line.a), lay.px(line.b));
                svg.open("line")
                    .num("x1", x1)
                    .num("y1", y1)
                    .num("x2", x2)
                    .num("y2", y2)
                    .empty();
            }
            for line in item.lines.iter().filter(|l| l.arrow) {
                self.arrowhead(svg, lay, line, item.style);
            }
            svg.close("g");
        }
        for item in &self.items {
            if item.dots.is_empty() {
                continue;
            }
            let colour = self.colour(item.style.colour);
            svg.open("g").attr("class", item.dot_class);
            match item.dots.first() {
                Some(dot) if dot.hollow => svg
                    .attr("fill", &surface)
                    .attr("stroke", &colour)
                    .num("stroke-width", 1.5),
                _ => svg
                    .attr("fill", &colour)
                    .attr("stroke", &surface)
                    .num("stroke-width", 1.5),
            };
            svg.children();
            Self::tooltip(svg, &item.label);
            for dot in &item.dots {
                let (cx, cy) = lay.px(dot.p);
                svg.open("circle")
                    .num("cx", cx)
                    .num("cy", cy)
                    .num("r", dot.radius)
                    .empty();
            }
            svg.close("g");
        }
        self.write_labels(svg, lay, &surface, &ink, &secondary);

        // Axes: tick marks, tick labels and titles.
        svg.open("g")
            .attr("class", "axis x")
            .attr("stroke", &muted)
            .num("stroke-width", 1.0)
            .attr("fill", &secondary)
            .num("font-size", TICK_SIZE)
            .attr("text-anchor", "middle")
            .children();
        for k in f.x_ticks.indices() {
            let x = x_at(k);
            svg.open("line")
                .num("x1", x)
                .num("y1", bottom)
                .num("x2", x)
                .num("y2", bottom + TICK_LENGTH)
                .empty();
            svg.open("text")
                .num("x", x)
                .num("y", bottom + X_LABEL_DROP)
                .attr("stroke", "none")
                .content()
                .text(&f.x_ticks.label(k))
                .close("text");
        }
        if !self.x_title.is_empty() || f.x_ticks.exponent != 0 {
            svg.open("text")
                .attr("class", "axis-title")
                .num("x", 0.5 * (left + right))
                .num("y", bottom + X_LABEL_DROP + AXIS_SIZE + 7.0)
                .num("font-size", AXIS_SIZE)
                .attr("stroke", "none");
            Self::axis_title(svg, &self.x_title, f.x_ticks.exponent);
        }
        svg.close("g");
        svg.open("g")
            .attr("class", "axis t")
            .attr("stroke", &muted)
            .num("stroke-width", 1.0)
            .attr("fill", &secondary)
            .num("font-size", TICK_SIZE)
            .attr("text-anchor", "end")
            .children();
        for k in f.t_ticks.indices() {
            let y = y_at(k);
            svg.open("line")
                .num("x1", left - TICK_LENGTH)
                .num("y1", y)
                .num("x2", left)
                .num("y2", y)
                .empty();
            svg.open("text")
                .num("x", left - TICK_LENGTH - 3.0)
                .num("y", y + 4.0)
                .attr("stroke", "none")
                .content()
                .text(&f.t_ticks.label(k))
                .close("text");
        }
        if !self.t_title.is_empty() || f.t_ticks.exponent != 0 {
            let (x, y) = (PAD + AXIS_SIZE - 2.0, 0.5 * (top + bottom));
            svg.open("text")
                .attr("class", "axis-title")
                .num("x", x)
                .num("y", y)
                .attr(
                    "transform",
                    &format!("rotate(-90 {} {})", fixed(x, 2), fixed(y, 2)),
                )
                .num("font-size", AXIS_SIZE)
                .attr("text-anchor", "middle")
                .attr("stroke", "none");
            Self::axis_title(svg, &self.t_title, f.t_ticks.exponent);
        }
        svg.close("g");

        if !self.title.is_empty() {
            svg.open("text")
                .attr("class", "title")
                .num("x", PAD)
                .num("y", f.title_base)
                .num("font-size", TITLE_SIZE)
                .attr("font-weight", "600")
                .attr("fill", &ink)
                .content()
                .text(&self.title)
                .close("text");
        }
        if !self.subtitle.is_empty() {
            svg.open("text")
                .attr("class", "subtitle")
                .num("x", PAD)
                .num("y", f.subtitle_base)
                .num("font-size", SUBTITLE_SIZE)
                .attr("fill", &secondary)
                .content()
                .text(&self.subtitle)
                .close("text");
        }

        if !f.legend.is_empty() {
            svg.open("g")
                .attr("class", "legend")
                .num("font-size", LEGEND_SIZE)
                .attr("fill", &ink)
                .children();
            let x = right + LEGEND_GAP;
            for (row, entry) in f.legend.iter().enumerate() {
                let y = top + 4.0 + (row as f64 + 0.5) * LEGEND_ROW;
                let colour = self.colour(entry.style.colour);
                if entry.swatch == Swatch::Dot {
                    svg.open("circle")
                        .num("cx", x + 0.5 * LEGEND_SWATCH)
                        .num("cy", y)
                        .num("r", TICK_RADIUS)
                        .attr("fill", &colour)
                        .attr("stroke", &surface)
                        .num("stroke-width", 1.5)
                        .empty();
                } else {
                    svg.open("line")
                        .num("x1", x)
                        .num("y1", y)
                        .num("x2", x + LEGEND_SWATCH)
                        .num("y2", y)
                        .attr("stroke", &colour)
                        .num("stroke-width", entry.style.width);
                    if let Some(pattern) = entry.style.dash.pattern() {
                        svg.attr("stroke-dasharray", pattern);
                    }
                    svg.empty();
                }
                svg.open("text")
                    .num("x", x + LEGEND_SWATCH + 8.0)
                    .num("y", y + 4.0)
                    .content()
                    .text(&entry.label)
                    .close("text");
            }
            svg.close("g");
        }
        svg.close("svg");
    }

    /// The text beside the dots. Names go first and are always written. A clock reading is left out where
    /// it would run into a dot or into text already placed, or off the plot, so crowded ticks lose numbers,
    /// not legibility.
    fn write_labels(&self, svg: &mut Svg, lay: &Layout, surface: &str, ink: &str, secondary: &str) {
        let dots = || self.items.iter().flat_map(|i| &i.dots);
        if dots().all(|d| d.text.is_empty()) {
            return;
        }
        let plot = [
            lay.plot_left,
            lay.plot_top,
            lay.plot_left + lay.plot_width,
            lay.plot_top + lay.plot_height,
        ];
        // Boxes as [left, top, right, bottom], in pixels.
        let overlap =
            |a: &[f64; 4], b: &[f64; 4]| a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3];
        let mut taken: Vec<[f64; 4]> = dots()
            .map(|d| {
                let (x, y) = lay.px(d.p);
                [x - d.radius, y - d.radius, x + d.radius, y + d.radius]
            })
            .collect();
        svg.open("g")
            .attr("class", "labels")
            .num("font-size", MARK_SIZE)
            .children();
        for quiet in [false, true] {
            for dot in dots().filter(|d| d.quiet == quiet && !d.text.is_empty()) {
                let (px, py) = lay.px(dot.p);
                let reach = dot.radius + 5.0;
                let width = text_width(&dot.text, MARK_SIZE);
                // Where the text goes for a placement: its anchor point, its anchor and its box.
                let place = |placement: Placement| {
                    let (x, y, anchor, left) = match placement {
                        Placement::Left => (px - reach, py + 4.0, "end", px - reach - width),
                        Placement::Above => (px, py - reach, "middle", px - 0.5 * width),
                        Placement::Below => (px, py + reach + 8.0, "middle", px - 0.5 * width),
                        _ => (px + reach, py + 4.0, "start", px + reach),
                    };
                    let area = [left, y - 0.9 * MARK_SIZE, left + width, y + 0.3 * MARK_SIZE];
                    (x, y, anchor, area)
                };
                let free = |area: &[f64; 4]| !taken.iter().any(|other| overlap(area, other));
                let inside = |area: &[f64; 4]| {
                    area[0] >= plot[0]
                        && area[1] >= plot[1]
                        && area[2] <= plot[2]
                        && area[3] <= plot[3]
                };
                let placement = match dot.placement {
                    // The first side with room inside the plot; failing that, the right.
                    Placement::Auto => [
                        Placement::Right,
                        Placement::Left,
                        Placement::Above,
                        Placement::Below,
                    ]
                    .into_iter()
                    .find(|&side| {
                        let area = place(side).3;
                        inside(&area) && free(&area)
                    })
                    .unwrap_or(Placement::Right),
                    other => other,
                };
                let (x, y, anchor, area) = place(placement);
                if quiet && !(free(&area) && inside(&area)) {
                    continue;
                }
                taken.push(area);
                // The text twice: first as a halo in the surface colour, then the letters.
                svg.open("text")
                    .num("x", x)
                    .num("y", y)
                    .attr("text-anchor", anchor)
                    .attr("fill", surface)
                    .attr("stroke", surface)
                    .num("stroke-width", 3.0)
                    .attr("stroke-linejoin", "round")
                    .content()
                    .text(&dot.text)
                    .close("text");
                svg.open("text")
                    .num("x", x)
                    .num("y", y)
                    .attr("text-anchor", anchor)
                    .attr("fill", if quiet { secondary } else { ink })
                    .content()
                    .text(&dot.text)
                    .close("text");
            }
        }
        svg.close("g");
    }

    /// The text of an axis title: the caller's words, then the power of ten taken out of the tick labels.
    /// Finishes the `<text>` element the caller opened.
    fn axis_title(svg: &mut Svg, title: &str, exponent: i32) {
        svg.content().text(title);
        if exponent != 0 {
            let power = if exponent < 0 {
                format!("\u{2212}{}", exponent.unsigned_abs())
            } else {
                exponent.to_string()
            };
            svg.text(if title.is_empty() { "×10" } else { " ×10" });
            svg.open("tspan")
                .num("dy", -4.0)
                .num("font-size", 9.0)
                .content()
                .text(&power)
                .close_inline("tspan");
        }
        svg.close("text");
    }

    /// A filled arrowhead at the middle of a line, pointing from its start to its end. Left out on a line too
    /// short to carry one.
    fn arrowhead(&self, svg: &mut Svg, lay: &Layout, line: &Line, style: Style) {
        let ((x1, y1), (x2, y2)) = (lay.px(line.a), lay.px(line.b));
        let (dx, dy) = (x2 - x1, y2 - y1);
        let length = libm::hypot(dx, dy);
        let size = 6.0 + 1.5 * style.width;
        if !(length >= 2.5 * size) {
            return;
        }
        let (ux, uy) = (dx / length, dy / length);
        let (mx, my) = (0.5 * (x1 + x2), 0.5 * (y1 + y2));
        let (tip, back, half) = (0.5 * size, -0.5 * size, 0.4 * size);
        svg.open("polygon")
            .attr("class", "arrow")
            .points([
                (mx + tip * ux, my + tip * uy),
                (mx + back * ux - half * uy, my + back * uy + half * ux),
                (mx + back * ux + half * uy, my + back * uy - half * ux),
            ])
            .attr("fill", &self.colour(style.colour))
            .attr("stroke", "none")
            .empty();
    }

    /// A sentence or two for readers who cannot see the picture.
    fn description(&self, f: &Furniture) -> String {
        let name = |title: &str, fallback: &str| {
            if title.is_empty() {
                fallback.to_string()
            } else {
                title.to_string()
            }
        };
        let mut text = format!(
            "Spacetime diagram: {} across, {} up.",
            name(&self.x_title, "position"),
            name(&self.t_title, "time")
        );
        if !self.subtitle.is_empty() {
            text.push(' ');
            text.push_str(&self.subtitle);
            if !self.subtitle.ends_with('.') {
                text.push('.');
            }
        }
        let mut names: Vec<&str> = f.legend.iter().map(|e| e.label.as_str()).collect();
        for dot in self.items.iter().flat_map(|i| &i.dots) {
            if !dot.quiet && !dot.text.is_empty() && !names.contains(&dot.text.as_str()) {
                names.push(&dot.text);
            }
        }
        if !names.is_empty() {
            text.push_str(" Shown: ");
            text.push_str(&names.join("; "));
            text.push('.');
        }
        text
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pixels_grow_right_and_down() {
        let layout = Layout {
            width: 300.0,
            height: 300.0,
            plot_left: 50.0,
            plot_top: 20.0,
            plot_width: 200.0,
            plot_height: 100.0,
            x: (-1.0, 3.0),
            t: (0.0, 2.0),
        };
        assert_eq!(layout.to_pixels(-1.0, 2.0), (50.0, 20.0));
        assert_eq!(layout.to_pixels(3.0, 0.0), (250.0, 120.0));
        assert_eq!(layout.to_pixels(1.0, 1.0), (150.0, 70.0));
    }
}
