# SRKit

Special relativity in 3+1 Minkowski spacetime from closed forms: piecewise-analytic worldlines, exact proper time, light-cone queries, an exact integer clock, flat ΛCDM seen from inside, and a small set of general-relativity scalars. Rust, with JavaScript (WebAssembly) and Python bindings.

It is the engine for *The Special Relativity Game* and is written to stand alone. The plan it follows is [`docs/engine/spacetime-package-plan.md`](../docs/engine/spacetime-package-plan.md). It lives in this repository for now and is laid out to move to its own (R71.7): nothing under `srkit/` depends on the game.

**Status: first complete version, 0.2.0, not published.** The licence (R24.6), the name (R24.8) and the public repository (R24.9) are open questions, so every crate is `publish = false` and there is no licence file. Code and text are AI-generated (Claude, building on a first core by Codex); see the generated-content note at the end.

## The crates

| Crate | What it is | `no_std` |
|---|---|---|
| [`srkit-core`](crates/srkit-core) | Segments, worldlines, light cones, itineraries, rides, optics, rockets, horizons, charts, the exact clock and exact causal predicates | yes; `alloc` only for owned worldlines |
| [`srkit-flrw`](crates/srkit-flrw) | Flat ΛCDM in closed form: scale factor, conformal time, the event horizon, reach, last contact, Hubble drag, the de Sitter limit | yes |
| [`srkit-gr`](crates/srkit-gr) | Closed-form scalars: Schwarzschild and Kerr time dilation and ISCO, Shapiro delay, deflection, Hawking radiation | yes |
| [`srkit-log`](crates/srkit-log) | An append-only event log read through light cones: each receiver gets events in the order their light arrives | no |
| [`srkit-edu`](crates/srkit-edu) | Minkowski diagrams as deterministic SVG | no |
| [`srkit-wasm`](crates/srkit-wasm) | The WebAssembly binding: a plain C ABI. [`js/srkit.src.js`](js/srkit.src.js) is its loader, [`js/srkit.d.ts`](js/srkit.d.ts) its types | no |
| [`srkit-py`](crates/srkit-py) | The Python binding (PyO3, abi3). Its own small workspace: it needs a Python toolchain | no |

## What the core does

**Units.** Geometry is in natural units with `c = 1`: choose a unit of time and measure lengths in the distance light covers in it (seconds and light-seconds; Julian years and light-years in the game). Velocities are fractions of `c`, rapidities are dimensionless, proper accelerations are rates of change of rapidity. `srkit_core::units` has the exact SI constants and `Units` converts both ways. Rapidity, not velocity, is a segment's state, so a ride at γ 10⁶ keeps every digit.

**Segments** ([`segment.rs`](crates/srkit-core/src/segment.rs)). Each gives position, velocity and proper time as explicit functions of coordinate time, and the inverse. Nothing is stepped.

| Kind | Motion | Proper time |
|---|---|---|
| `Rest`, `Inertial` | Still, or at constant velocity | Linear |
| `Hyperbolic` | Constant proper acceleration along the motion | Closed form (`log1p`) |
| `Steered` | Constant proper acceleration at any angle, the direction fixed aboard | Closed form |
| `Burn` | Constant thrust while the mass falls, for any exhaust speed | Closed form by proper time; a bracketed root by coordinate time |
| `Circular` | A circle at constant speed | Linear |
| `Null` | A light-like jump | None passes |

**Worldlines and light cones** ([`worldline.rs`](crates/srkit-core/src/worldline.rs)). `send(event)` is the moment on a worldline whose light reaches an observer: what that observer sees. `recv(event)` is when light from an event reaches the worldline. Both are the unique root of an increasing function: a closed form on rest, inertial and null segments, a bracketed Newton iteration (about five evaluations) on the others. A ship that accelerates for ever has a horizon, and the search asks the horizon's closed form first, so "never" is an answer and not a timeout.

**Itineraries** ([`plan.rs`](crates/srkit-core/src/plan.rs)). Runs that join with no jump in velocity: a straight line between two rapidities with a speed cap and separate boost and brake accelerations; a brake to rest; entering a circular orbit along a tangent; leaving one. [`trip.rs`](crates/srkit-core/src/trip.rs) is the Navigator's ride in null coordinates (legs, a γ cap, stopovers, a return), its inverse and its light-cone intersection with any event. [`reach.rs`](crates/srkit-core/src/reach.rs) is where a ship can come to rest on a proper-time budget.

**Seeing** ([`optics.rs`](crates/srkit-core/src/optics.rs), [`kinematics.rs`](crates/srkit-core/src/kinematics.rs), [`frames.rs`](crates/srkit-core/src/frames.rs)). Aberration and Doppler in a ship's axes, the rate a distant clock appears to run, boosts by a rapidity vector, radar and light-cone coordinates of an event for any observer, diffraction and link budgets. [`rindler.rs`](crates/srkit-core/src/rindler.rs): the horizon of an accelerating ship and the last moment a signal from a given place can still catch it.

**Exact time** ([`time.rs`](crates/srkit-core/src/time.rs), [`causal.rs`](crates/srkit-core/src/causal.rs)). Stored time is an unsigned count of milliseconds in 384 or 512 bits, with checked arithmetic and integer modular phase. Causal decisions are made in integers: with positions rounded to a one-metre grid, "has this signal arrived by tick n?" is `(c·Δt)² ≥ 10⁶·|Δx|²`, which holds or fails with no rounding. A moving receiver's arrival is found in floating point and then settled to the first tick that passes the exact check. This is the "exact predicates, inexact constructions" rule of computational geometry (plan §3.1).

## Build and test

Rust is pinned by `rust-toolchain.toml` (1.85.1, with rustfmt, Clippy and the wasm32 target) and dependencies by `Cargo.lock`. From this directory:

```sh
cargo test --workspace --locked                 # all crates, debug
cargo test --workspace --release --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo fmt --all -- --check
cargo run --example trip                        # a 1 g ride to Proxima, start to finish
cargo run --release --example clock_widths      # U384 against U512
```

The JavaScript build (from the repository root; Python 3 and Node):

```sh
python3 -I srkit/tools/build_js.py              # srkit/dist/srkit.js, srkit.wasm, srkit.d.ts
node --expose-gc srkit/tools/test_js.js         # the binding's own tests
node srkit/tools/bits.js                        # native and wasm32 give the same bits
node srkit/tools/bench_clocks.js                # U384 against U512 in WebAssembly
```

One script runs all of the above, and with `--apps` the two applications' checks as well:

```sh
srkit/tools/check_all.sh --apps
```

`build_js.py` remaps source paths, so one commit gives the same module bytes on any machine; the module's SHA-256 is in `SRKit.info.sha256` and is shown by both applications. The applications carry their own copy of the built script (see below).

The Python build is in [`crates/srkit-py/README.md`](crates/srkit-py/README.md).

## Using it

Rust:

```rust
use srkit_core::{kinematics::Event, plan::{self, Start}, segment::Segment, units::G0_PER_YEAR, worldline::WorldlineBuf, Vec3};

// years and light-years: a 1 g ride to Proxima
let proxima = Vec3::new(4.2465, 0.0, 0.0);
let ride = plan::line(Start { t: 2100.0, tau: 0.0, x: Vec3::ZERO }, proxima, G0_PER_YEAR, f64::INFINITY, 0.0, 0.0, G0_PER_YEAR, None).unwrap();
assert!((ride.end.tau - 3.542).abs() < 1e-3);          // years aboard

let mut segments = vec![Segment::rest(2000.0, 2100.0, Vec3::ZERO)];
segments.extend(ride.pieces.iter().map(|p| p.seg));
segments.push(Segment::rest(ride.end.t, f64::INFINITY, proxima));
let ship = WorldlineBuf::new(segments, 2100.0, 0.0).unwrap();

// what Earth sees of the ship in 2105: a query on Earth's past light cone
let (t_emitted, seen) = ship.view().seen_from(Event::new(2105.0, Vec3::ZERO)).unwrap();
```

JavaScript (one script, no build step in the page):

```js
const g = SRKit.units.G_PER_YEAR;                                   // 1 g per year
const ride = SRKit.plan.line({ t: 2100, tau: 0, X: [0, 0, 0] }, [4.2465, 0, 0], g);
const ship = SRKit.wl.make([SRKit.seg.mk.rest(2000, 2100, [0, 0, 0]), ...ride.pieces.map(p => p.g),
                            SRKit.seg.mk.rest(ride.end.t, Infinity, ride.end.X)], 2100, 0);
const tSeen = SRKit.wl.send(ship, 2105, [0, 0, 0]);                 // when the light Earth sees in 2105 left the ship
SRKit.time.arrived('3' + '0'.repeat(110), '3' + '0'.repeat(98) + '134009348401', [4.0175e16, 0, 0]);   // exact, at any epoch
```

## How it is checked

197 Rust tests (core 41, cosmology 41, gravity 50, event log 8, diagrams 54, binding 3), 951 checks of the JavaScript binding and 154 of the Python binding. What they are:

- **Golden vectors.** `tools/golden_core.py` computes 610 cases from their textbook formulas (or a quadrature, or a root found by mpmath) at 60 significant digits, sharing no code with the crate. `tests/golden.rs` holds each family to a stated relative tolerance. Measured worst errors: 1e-15 to 5e-15 for segments, rides, optics, boosts, burns and Kepler's equation; 1e-11 where the input is a velocity near c; 6e-12 for reach (its iteration stops at 1e-12 by design).
- **The physics oracle.** `tools/physics-checks/oracle.csv` is the project's table of reference numbers. `srkit-core` recomputes 47 of its rows through the public API, `srkit-flrw` and `srkit-gr` the rows that are theirs. Rows that are game design are left to `oracle_check.py`, which still reproduces every row.
- **Properties** (proptest): proper time inverts; nothing outruns light; velocity is the derivative of position; light-cone roots lie on the cone; the interval is invariant under boosts; steered thrust holds its proper acceleration and agrees with the straight form when aligned; burns obey the rocket equation; itineraries join in position and velocity and end on their target; the event log gives one history however it is read.
- **Exact causality.** On the oracle's adversarial test (arrivals a hair after a tick, DT-14), `ceil()` of the floating-point answer is a tick early 2,678 times in 100,000; the settled tick is never early and never late.
- **The applications.** The Navigator and alpha 3.3e3 run on it, and each has a parity test against the JavaScript it replaced: [`docs/verification/srkit-ports.md`](../docs/verification/srkit-ports.md).

## Determinism

Every transcendental goes through the pure-Rust `libm`; there is no `unsafe` outside the binding crates' export attributes, no threads, no hash-map iteration and no fused multiply-add. One build therefore gives the same bits on any IEEE 754 platform, and the WebAssembly build gives the same bits in every browser, which JavaScript's own `Math` does not promise.

`tools/bits.js` checks this where it can: it runs a fixed list of calls natively and through the WebAssembly build and compares bit patterns. On x86-64 the two agree on every result. `srkit-flrw` pins 40 of its results to the bit in a test.

Not yet done from the plan's contract (§3): the same comparison on aarch64, a pinned golden file of bit patterns in CI for every crate, and an error bound in units in the last place stated for each function. The golden tests here compare to a tolerance. Until those are in place the contract is not frozen and a release is 0.x.

## Limits

- Geometry is binary64. An intersection is a construction, good to about 1e-15 of the scale of the problem; only [`causal`] decides exactly, and only on the grid.
- By coordinate time, a segment that brakes from rapidity φ₀ resolves its last moments to rounding × sinh φ₀. Evaluate by proper time near the turnaround, as the itineraries do.
- `Burn` is collinear: thrust along or against the motion.
- Circular segments use a floating-point phase. The clock module reduces an integer-millisecond period exactly; joining the two for deep-time orbits is not done.
- `srkit-flrw` is matter and Λ only. Radiation is out of scope for the closed forms.
- Kerr orbits off the equator, geodesics and imaging are out of scope for `srkit-gr` by design.
- The Python binding covers the core. `srkit-flrw`, `srkit-gr` and `srkit-log` are not bound in Python yet; the JavaScript binding has the core, cosmology and gravity, and not the log or the diagrams.

## Rough edges in the API

Found by writing three bindings against the core, and not yet smoothed:

- "No answer" comes back four ways: `Option` from a segment's light-cone root, `±inf` from a worldline's, NaN from an orbit's, `None` from `seen_from`.
- The plain segment constructors do not validate (evaluation is total, and applications build segments in hot loops). `Segment::validate` checks one from outside; `Steered` and `Burn` normalise their directions, the others trust them.
- A plan that starts with a boost, used as a worldline's first segment, continues backwards as a ship braking in from infinity. Put a `Rest` before it.
- `Trip` caps the Lorentz factor; `plan::line` caps the rapidity. `plan::brake` and `orbit_exit` return tuples where `line` and `orbit_entry` return structs.
- `Duration` has no ordering; `modulo_millis` takes a period of at most 64 bits; a `Timestamp` is built from decimal text or a `bnum` integer, so a binding needs `bnum` at the same version.

## Layout

```
srkit/
  Cargo.toml, Cargo.lock, rust-toolchain.toml
  crates/            the seven crates
  js/                srkit.src.js (the loader) and srkit.d.ts
  tools/             build_js.py, test_js.js, bits.js, bench_clocks.js, golden_core.py, check_all.sh
  dist/              built output (not tracked)
```

Applications that carry a built copy: `prototypes/alpha-3-3e3/shared/srkit.js` and, embedded in the page, `docs/design/mockups/relativistic-navigator.html`. Rebuild them with:

```sh
python3 -I srkit/tools/build_js.py --copy-to prototypes/alpha-3-3e3/shared/srkit.js
python3 -I docs/design/mockups/src/build_srkit.py
python3 -I docs/design/mockups/src/build_nav.py docs/design/mockups/src/navigator.src.html docs/design/mockups/src/stars.json docs/design/mockups/relativistic-navigator.html --standalone
```

## Generated content

Everything under `srkit/` is AI-generated: the Rust, the tests, the JavaScript and Python bindings, the tools and this file. The exact-clock module and several numerical forms come from the first core written by Codex on `codex/srkit-core`; the rest was written by Claude on `claude/busy-volta-gj61zs`. The record is [`docs/verification/branch-review-2026-11-16.md`](../docs/verification/branch-review-2026-11-16.md).
