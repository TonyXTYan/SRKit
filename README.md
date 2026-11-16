# SRKit: first Rust core

Experimental `no_std`, allocation-free special relativity primitives, built on
the draft [package plan](../docs/engine/spacetime-package-plan.md). This directory
keeps the physics separate from the game for later extraction into its own
repository. Publication, licensing, clock-width selection and the full P0
determinism contract remain open; the crate has `publish = false`.

## Build and check

Install Rust using rustup. `rust-toolchain.toml` pins Rust 1.85.1, rustfmt,
Clippy and the WebAssembly target. From this directory:

```sh
cargo test --workspace --locked
cargo test --workspace --release --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo fmt --all -- --check
cargo build --workspace --release --locked
cargo build --workspace --target wasm32-unknown-unknown --release --locked
cargo run --example trip --locked
```

In the prepared cloud environment, first set:

```sh
export CARGO_HOME=/workspace/.cargo
export RUSTUP_HOME=/workspace/.rustup
export PATH="$CARGO_HOME/bin:$PATH"
```

Tests read the existing `../tools/physics-checks/oracle.csv`, the single source
of oracle values. Keep that file with the checkout. Run the independent Python
oracle from the repository root using the scientific Python environment.

## Implemented

- SI units, full 3D vectors, Minkowski four-vectors, Lorentz boosts, Doppler
  factors and aberration. Photon directions always follow propagation.
- Bounded inertial, arbitrary-angle constant proper acceleration, circular and
  null segments, analytic proper-time evaluation and inversion. Acceleration
  acts along a direction fixed aboard, in the initial four-velocity/thrust
  plane. Worldlines validate position, time and velocity joins.
- Retarded source states and future light-cone intersections using bounded
  monotone root solving. Null continua choose the latest retarded or earliest
  future event. Queries outside a segment return an error.
- Symmetric flip-and-burn trips, photon-rocket mass ratios, circular gamma and
  a directly evaluated light-time excess that survives deep-time cancellation.
- Both U384 and U512 integer millisecond timestamps, signed durations, checked
  arithmetic, decimal interchange and integer modular phase.
- Exact causal round-up **for stationary integer-metre distances only**, with
  checked full-width arithmetic and tests rejecting one-tick-early arrivals.

Tests cover 33 of the 111 unique shared oracle rows. The other 78 are explicitly
outside the implemented oracle subset, including GR, cosmology, gameplay time
scales and the still-unimplemented high-rapidity segment path. Property tests
exercise boost invariance, analytic clocks, numerical light-cone brackets,
integer clock arithmetic and exact stationary arrivals.

## Numerical contract and remaining work

Exact integer timestamps do not make floating-point geometry exact. Subtract
epochs in integer arithmetic before using local coordinates. Local clock
conversion accepts at most ±2^53 ms; its division to seconds is rounded as
binary64. Segment events use local seconds/metres and reject non-finite values,
unrepresentable intervals and timelike velocities that round to light speed.
Very high rapidities are not yet supported as full segment states; trip scalar
formulas have a wider domain.

Circular segments evaluate a bounded local floating phase. The clock module can
reduce an *integer-millisecond period* exactly, but that is not an exact physical
representation of every circular period. These pieces are not yet a coupled
deep-time periodic segment implementation.

All transcendentals use pure-Rust `libm`; there is no unsafe code. This is not a
frozen cross-platform bitwise determinism guarantee. Golden vectors, execution
under wasm/aarch64, certified moving-worldline causal predicates, Rindler-horizon
queries, mass-burning segments, width benchmarks and all P1/P2 modules and
language bindings remain future work. A wasm library build verifies target
compatibility; it does not create a JavaScript binding or exercise browser code.

The review and validation record is in
[SRKit core review](../docs/verification/srkit-core-review.md).
