//! # srkit-edu
//!
//! Minkowski (spacetime) diagrams of [`srkit_core`] worldlines, written as SVG.
//!
//! A [`Diagram`] is a 1+1D picture: time up, one spatial axis across. Give it a [`View`] (the window, the size
//! of the plot, the axis in space and, if wanted, a moving [`Frame`] to draw in), add things to it, and call
//! [`Diagram::render`]:
//!
//! * worldlines, sampled so that the drawn line is within a quarter of a pixel of the true curve;
//! * events, light cones and light signals between two events;
//! * the lines of simultaneity and of constant position of any inertial frame (the tilted axes of a moving
//!   observer);
//! * ticks of proper time along a worldline;
//! * the Rindler horizon of an accelerating ship;
//! * the light ray by which an observer sees another body.
//!
//! [`twin_paradox`] and [`rindler_chase`] are two complete diagrams built this way.
//!
//! ```
//! use srkit_core::{kinematics::Event, segment::Segment, worldline::WorldlineBuf, Vec3};
//! use srkit_edu::{Cone, Diagram, Style, TickNumbers, View};
//!
//! // A ship leaving the origin at 0.6 c, in years and light-years.
//! let ship = Segment::inertial(0.0, 10.0, Vec3::ZERO, Vec3::new(0.6, 0.0, 0.0));
//! let ship = WorldlineBuf::new([ship], 0.0, 0.0).unwrap();
//!
//! let mut diagram = Diagram::new(View::equal_scales((-2.0, 8.0), (0.0, 10.0), 40.0))?;
//! diagram.title("A ship at 0.6 c").axis_titles("x (ly)", "t (yr)");
//! diagram
//!     .light_cone(Event::new(0.0, Vec3::ZERO), Cone::Future, "Light from the launch", Style::light())?
//!     .worldline(&ship.view(), "Ship", Style::series(0))?
//!     .proper_time_ticks(&ship.view(), 1.0, TickNumbers::Right, "Ship clock, yr", Style::series(0))?;
//! let svg = diagram.render();
//! assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
//! # Ok::<(), srkit_edu::Error>(())
//! ```
//!
//! ## Units and frames
//!
//! Natural units, as in `srkit-core`: `c = 1`, so one unit of time is one unit of length and, with equal
//! scales on the two axes ([`View::equal_scales`]), light runs at 45°. The diagram knows no unit names; put
//! them in the axis titles.
//!
//! Everything is handed over in the lab frame, the frame the worldlines are stored in. A view with a
//! [`Frame`] shows the same scene from an observer moving along the axis: every event is boosted with
//! [`srkit_core::kinematics::boost_rapidity`] about the frame's origin event before it is drawn.
//!
//! The diagram has one spatial axis. A body off the axis, or moving across it, is projected onto it; light
//! between two events that are not both on the axis is then drawn steeper than 45°.
//!
//! ## Output
//!
//! One self-contained SVG document: no scripts, no style sheets, no external fonts, images or links. Every
//! label is XML-escaped. Every coordinate is written with two decimals, tick labels are built from integers,
//! and all transcendental functions go through the pure-Rust `libm`, so the same calls give the same bytes on
//! every platform. There is no hash map and no clock anywhere in the crate.
//!
//! Two colour schemes, [`Theme::Light`] and [`Theme::Dark`]. The palette is fixed: see [`Theme`] for what was
//! checked. Colour is never the only cue, since series differ in dash pattern as well and every labelled line
//! is named in the legend.
#![forbid(unsafe_code)]
// Comparisons are written `!(x > 0.0)` on purpose: unlike `x <= 0.0` they also reject a NaN.
#![allow(clippy::neg_cmp_op_on_partial_ord, clippy::nonminimal_bool)]
#![warn(missing_docs)]

mod axis;
mod diagram;
mod error;
mod geometry;
mod sample;
mod scenes;
mod style;
mod svg;
mod view;

pub use diagram::{Cone, Diagram, Layout, Placement, TickNumbers, MAX_MARKS, NULL_TOLERANCE};
pub use error::Error;
pub use sample::{CURVE_TOLERANCE, MAX_POINTS};
pub use scenes::{rindler_chase, rindler_chase_diagram, twin_paradox, twin_paradox_diagram};
pub use style::{Colour, Dash, Style, Theme};
pub use view::{Frame, View, MAX_PIXELS};

/// Compiles and runs the example in the README as a doc test.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
