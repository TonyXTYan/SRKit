# srkit-gr

Closed-form general-relativity scalars for SRKit: numbers that have an exact formula, evaluated directly.
`no_std`, no allocation, no `unsafe`, every transcendental through the pure-Rust `libm`.

It is the `srkit-gr` layer of the [package plan](../../../docs/engine/spacetime-package-plan.md): "a small
set of closed-form scalars, each cross-checked".

## What it provides

All arguments and results are in SI units (kg, m, s, rad, J, W, K); ratios are dimensionless. Radii are
coordinate radii (Schwarzschild areal radius, Boyer-Lindquist radius). The spin of a Kerr hole is the
dimensionless `a* = Jc/(GM²)` in `[−1, 1]`, and orbits are `Orbit::Prograde` or `Orbit::Retrograde`
relative to the hole's rotation. Functions are total: outside the stated domain they return NaN.

| Module | Functions |
|---|---|
| crate root | `gravitational_parameter` (GM), `mass_from_gravitational_parameter`, `gravitational_radius` (GM/c²), `gravitational_time` (GM/c³) |
| `schwarzschild` | `schwarzschild_radius`, `photon_sphere_radius`, `marginally_bound_radius`, `isco_radius`; `static_dtau_dt` = √(1 − r_s/r); `gravitational_redshift` between two static radii; `circular_orbit_dt_dtau` = 1/√(1 − 3GM/(rc²)); `circular_orbit_angular_velocity` |
| `kerr` | `outer_horizon_radius`, `inner_horizon_radius`, `ergosphere_radius`, `horizon_angular_velocity`, `irreducible_mass`, `extractable_energy_fraction`; `static_dtau_dt`; `zamo_angular_velocity` (frame dragging) and `zamo_dtau_dt`; `photon_orbit_radius`, `marginally_bound_radius`, `isco_radius` (Bardeen, Press and Teukolsky 1972); `isco_specific_energy`, `isco_specific_angular_momentum`, `isco_dt_dtau`, `accretion_efficiency`; `circular_orbit_dt_dtau`, `circular_orbit_specific_energy`, `circular_orbit_specific_angular_momentum`, `circular_orbit_angular_velocity` at any radius |
| `weak_field` | `shapiro_delay` (two legs with an impact parameter), `shapiro_delay_leg`, `shapiro_delay_between` (two points and their separation); `deflection_angle` = 4GM/(bc²); `einstein_angle`; clock-rate offsets `static_rate_offset` = −GM/(rc²), `circular_orbit_rate_offset` = −3GM/(2rc²), `orbit_vs_surface_rate_offset` (the GPS sum) |
| `hawking` | `temperature`, `entropy`; `power`, `mass_loss_rate`, `lifetime`, `mass_for_lifetime`, each taking an `EmissionModel` |
| `newtonian` | `binding_energy_uniform_sphere`, `circular_speed`, `escape_speed`, `circular_period`, `synodic_period`, `hohmann_transfer_time`, `hohmann_delta_v` |

Every function's documentation gives its formula, the units of each argument and of the result, the domain
where it is valid, a literature reference and its accuracy.

### Hawking emission models

The temperature and entropy of a hole follow from its mass. Its power and lifetime also depend on which
particles it emits, so those functions take an explicit model, `P = α ħc⁶/(G²M²)`:

| `EmissionModel` | α | Meaning |
|---|---|---|
| `PhotonBlackbody` | 1/(15360π) = 2.072e-5 | Photons only, black body of the horizon area, no greybody factors. The textbook `t = 5120πG²M³/(ħc⁴)`. This is the model behind the oracle rows SR-64, SR-65, DT-21 and DT-22 |
| `PagePhotonsGravitons` | 3.740e-5 | Photons and gravitons with Page's (1976) greybody factors: the realistic model for any hole too cold to emit massive particles |
| `MacGibbon { f }` | 2.793e-4 × f | MacGibbon's (1991) `dM/dt = −5.34e16 f(M)/M²` kg/s; the caller supplies the species factor `f(M)` |
| `Coefficient(α)` | α | Anything else |

The crate does not compute which species a hole of a given mass emits (the `f(M)` of the realistic
emission decided in R25.5). It gives the scalar once the model is chosen. A lifetime holds α constant.

## What is out of scope

- **Geodesics and orbits that are not circular and equatorial**: integrating motion in a metric belongs to
  [EinsteinPy](https://github.com/einsteinpy/einsteinpy) (numerical) and to the Black Hole Perturbation
  Toolkit's [KerrGeodesics](https://github.com/BlackHolePerturbationToolkit/KerrGeodesics) and
  [KerrGeoPy](https://github.com/BlackHolePerturbationToolkit/KerrGeoPy) (analytic Kerr orbits).
- **Imaging and ray tracing** (shadows, lensed images, strong deflection): Gyoto, ipole, kgeo, AART.
- **Metric and tensor algebra**: EinsteinPy's symbolic module.
- **Cosmology**: `srkit-flrw`. **Special relativity**, including the kinetic energy of oracle row SR-63:
  `srkit-core`.
- Charged holes, perturbations and gravitational waves, extended or oblate bodies, post-Newtonian orders
  beyond the first, and the particle physics of Hawking emission.

## How it was checked

Run from the workspace directory `srkit/`:

```sh
cargo test -p srkit-gr
cargo clippy -p srkit-gr --all-targets --no-deps -- -D warnings
cargo fmt -p srkit-gr -- --check
cargo build -p srkit-gr --target wasm32-unknown-unknown --release
```

- **Physics oracle** (`tests/oracle.rs`): the 24 gravity rows of `tools/physics-checks/oracle.csv` (SR-40 to
  SR-43, SR-47 to SR-49, SR-53 to SR-65, DT-19 to DT-22), read from the CSV and compared at each row's own
  tolerance. Nine rows are quoted to two figures and lie more than 1% from the true value; for exactly
  those the test accepts a value that rounds to the stated one, as `oracle_check.py` does (its ROUNDED
  status). It also reproduces the numbers `hawking_scaling.py` prints and the alternative emission models
  quoted in the oracle's notes.
- **Independent references** (`tools/crosscheck.py`, `tests/crosscheck.rs`): 1,400 values computed with
  mpmath at 120 digits, by a different route from the crate wherever one exists (metric components and root
  finding for Kerr, numerical integrals for the Shapiro delay, the deflection, the binding energy and the
  Hawking lifetime, `σAT⁴` for the photon black body, the vis-viva equation for Hohmann burns). Every public
  function is covered. The tolerance is `2e-15` except where the quantity is ill conditioned in its
  arguments, and the worst error at that tolerance is `7.4e-16`. Regenerate with
  `python crates/srkit-gr/tools/crosscheck.py`; `--check` fails if the checked-in CSV is stale. The script
  pins the constants of `srkit_core::units` and the test compares them bit for bit.
- **Limits and textbook values** (`tests/limits.rs`): Kerr reduces to Schwarzschild at zero spin; the ISCO
  is at 1 GM/c² (prograde) and 9 GM/c² (retrograde) for the extremal hole; efficiencies 1 − √(8/9) and
  1 − 1/√3; the Sun-grazing Earth-Mars Shapiro delay (247 µs round trip); 1.75″ deflection at the solar limb;
  the GPS clock offset (+38.5 µs/day for a point-mass Earth, +38.6 with the geoid potential); the first law
  `dE = T dS`; Earth-to-Mars Hohmann numbers; NaN outside every domain.
- **Properties** (`tests/properties.rs`, proptest): the ISCO radius is monotonic in spin, efficiencies and
  clock rates stay within their bounds, the characteristic radii are ordered, redshifts compose, the two
  forms of the Shapiro delay agree, Hawking quantities scale as powers of the mass, and no function panics
  on any binary64 input.

- **Same bits natively and in WebAssembly** (by hand, 16 Nov 2026, not automated): the 1,400 reference inputs
  were run through a release build on x86-64 and through the `wasm32-unknown-unknown` build under Node 22
  (V8); every result was bit-identical.

Not yet done: the cross-checks against EinsteinPy and KerrGeoPy in CI that the plan calls for (the references
here are mpmath only), and golden vectors run automatically on aarch64 and wasm32.

## Known limits

- A spin is a binary64 number, so `1 − a*` is resolved only to `1.1e-16`. Quantities on the prograde ISCO of
  a nearly extremal hole grow as `(1 − a*)^(−1/3)` and inherit that: `1 − 1e-14` is stored as
  `1 − 9.992e-15`, which moves `dt/dτ` by 0.03%. A caller that needs holes closer to extremal than about
  `1 − 1e-12` should have an interface in terms of `1 − a*`, which this crate does not provide.
- Clock rates very close to a horizon are as well determined as `r − r₊`; see the module documentation.
- Weak-field functions are first order in `GM/(rc²)` for a point mass.
