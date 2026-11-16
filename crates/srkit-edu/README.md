# srkit-edu

Minkowski (spacetime) diagrams of `srkit-core` worldlines, written as SVG. Time runs up the page and one
spatial axis runs across. The crate depends on `srkit-core` and the workspace's `libm`, and on nothing else.

This is the SVG export of the `srkit-edu` row in the
[package plan](../../../docs/engine/spacetime-package-plan.md). The notebooks and the browser playground of
that row are not here.

## What it draws

- **Worldlines**, from a `srkit_core::worldline::Worldline`. Straight segments (rest, inertial, light-like)
  take two points. Curved ones (constant proper acceleration along or across the motion, a burn that spends
  mass, a circle) are cut until the drawn line is within a quarter of a pixel of the true curve. Nothing is
  evaluated before the worldline's first segment starts or after its last one ends.
- **Events**, as dots with a label.
- **Light cones** of an event, past, future or both.
- **Light signals** between two events, with an arrowhead. The two events must be joined by light, or the
  call returns an error.
- **Another frame's lines**: lines of simultaneity and of constant position of any inertial frame moving
  along the axis, its two axes, or its whole coordinate grid. These are the tilted axes of a moving observer.
- **Proper-time ticks** along a worldline, every Δτ, with the clock reading beside each if wanted.
- **The Rindler horizon** of a segment of constant proper acceleration.
- **What is seen**: the light ray from a source's worldline to an observer, found with `Worldline::send`.
- Axes with ticks at steps of 1, 2 or 5 times a power of ten, the caller's axis titles (the diagram knows no
  units), a title, a subtitle and a legend beside the plot.

A `View` sets the window, the size of the plot in pixels, the axis in space, and optionally a `Frame`: the
same scene is then drawn as an observer moving along the axis finds it, with every event boosted about the
frame's origin event. Light is at 45° in every frame when the two axes share a scale.

`twin_paradox(distance, acceleration)` and `rindler_chase(acceleration)` return two complete diagrams.

## Use

```rust
use srkit_core::{kinematics::Event, segment::Segment, worldline::WorldlineBuf, Vec3};
use srkit_edu::{Cone, Diagram, Error, Style, TickNumbers, View};

fn main() -> Result<(), Error> {
    // A ship leaving the origin at 0.6 c, in years and light-years.
    let ship = Segment::inertial(0.0, 10.0, Vec3::ZERO, Vec3::new(0.6, 0.0, 0.0));
    let ship = WorldlineBuf::new([ship], 0.0, 0.0).map_err(|_| Error::NonFinite)?;
    // 40 pixels to the year and to the light-year, so light runs at 45 degrees.
    let mut diagram = Diagram::new(View::equal_scales((-2.0, 8.0), (0.0, 10.0), 40.0))?;
    diagram.title("A ship at 0.6 c").axis_titles("x (ly)", "t (yr)");
    diagram
        .light_cone(Event::new(0.0, Vec3::ZERO), Cone::Future, "Light from the launch", Style::light())?
        .worldline(&ship.view(), "Ship", Style::series(0))?
        .proper_time_ticks(&ship.view(), 1.0, TickNumbers::Right, "Ship clock, every year", Style::series(0))?;
    print!("{}", diagram.render());
    Ok(())
}
```

Three examples print SVG to standard output:

```sh
cargo run -p srkit-edu --example twin_paradox > twins.svg       # Proxima Centauri and back at one gravity
cargo run -p srkit-edu --example rindler_chase -- 0.5 dark > chase.svg
cargo run -p srkit-edu --example frames -- ship > ship.svg      # one scene, drawn in the ship's frame
```

## Output

One self-contained SVG document: no scripts, no style sheets, no external fonts, images or links. Every
label is XML-escaped. Each line carries its label as an SVG `<title>`, which browsers show on hover.

**The output is deterministic, byte for byte.** The same calls in the same order give the same string on
every platform: coordinates are written with two decimals, tick labels are built from integers, every
transcendental function goes through the pure-Rust `libm`, and the crate uses no hash map, clock or random
number. `tests/golden/twin_paradox.svg` is the checked-in rendering of `twin_paradox(4.2465, 1.0323)`; a test
compares it with the output byte for byte.

No function panics on its input. Bad input is an `Error`: an empty view, a number that is not finite, a
signal that is not light-like, a step that would draw thousands of ticks.

## Colour

Two themes, `Theme::Light` and `Theme::Dark`. The palette is fixed and was checked by computation: in both
themes the four series colours and the colour for light have a contrast of at least 3:1 against the surface,
every pair of them stays apart under simulated protanopia and deuteranopia, and text is at least 7:1. Colour
is never the only cue. `Style::series(n)` gives each series its own dash pattern as well, light is drawn
dashed, and every labelled line is named in the legend, so a diagram can be read in greyscale.

## Limits

- One spatial axis. A body off the axis is projected onto it, and light between events that are not both on
  the axis is then drawn steeper than 45°.
- Text is not measured, because the output names only system fonts. Margins and the legend are laid out from
  an estimate of the width of each label that leans wide.
- Labels beside dots avoid each other in a simple way: a clock reading that would run into other text is left
  out. Lines are not moved out of the way of labels.

## Checks

```sh
cargo test -p srkit-edu
cargo clippy -p srkit-edu --all-targets -- -D warnings
cargo fmt -p srkit-edu -- --check
```

To replace the golden file after a change that is meant, look at the new picture first, then:

```sh
cargo run -p srkit-edu --example twin_paradox > crates/srkit-edu/tests/golden/twin_paradox.svg
```
