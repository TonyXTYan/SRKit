# srkit-flrw

Flat ΛCDM cosmology in closed form, as seen by an observer inside the universe. `no_std`, no allocation, no
`unsafe`, every transcendental through the pure-Rust `libm`.

It answers the questions a game set in an expanding universe asks: how far can a signal sent now ever get, when
is the last moment to call a distant galaxy, when does light sent across the void arrive, how fast can a rocket
go against Hubble drag, and how cold does the sky get.

```rust
use srkit_flrw::Cosmology;

let u = Cosmology::PLANCK_2018;      // Julian years and light-years
let now = u.age();                   // 13.80e9 yr
let horizon = u.event_horizon(now);  // 16.68e9 ly, comoving

// The Virgo cluster, 53.8 million light-years away (comoving):
let virgo = 5.38158e7;
let arrival = u.light_arrival_time(now, virgo); // Some(now + 53.9e6 yr)
let last = u.last_contact_time(virgo);          // Some(now + 100.3e9 yr)
let lost = u.light_arrival_time(now, 20e9);     // None: beyond the event horizon
```

## What it provides

All functions are methods of `Cosmology`, a flat universe of matter and a cosmological constant fixed by `H0`
and `Ω_m` (`Ω_Λ = 1 − Ω_m`). Presets: `Cosmology::PLANCK_2018` (`H0` = 67.4 km/s/Mpc, `Ω_m` = 0.315; Julian
years and light-years) and `Cosmology::PLANCK_2018_SI` (the same, in seconds and light-seconds). Constructors:
`per_year`, `per_second`, `with_time_unit`, `from_km_s_mpc`, `from_h_lambda_per_year`.

| Group | Functions |
|---|---|
| Background | `scale_factor`, `ln_scale_factor`, `time_at_scale_factor`, `time_at_ln_scale_factor`, `hubble`, `hubble_at_scale_factor`, `age`, `time_at_redshift`, `lookback_time`, `redshift_between` |
| Conformal time | `conformal_time`, `conformal_time_at_scale_factor`, `conformal_time_limit`, `conformal_time_remaining`, `time_at_conformal_time`, `time_at_conformal_time_remaining` |
| Distances and horizons | `comoving_distance`, `comoving_distance_at_redshift`, `particle_horizon`, `particle_horizon_proper`, `event_horizon`, `event_horizon_proper` |
| Reach | `reach`, `reach_volume_ratio`, `reachable_share_of_observable`, `reach_mass` |
| Last contact | `last_contact_time`, `last_light_time`, `last_reply_time` |
| Light cones | `light_reception` (conformal chart), `light_arrival_time`, `light_emission_time` |
| Hubble drag | `coasting_proper_velocity`, `coasting_rapidity`, `coasting_gamma`, `scale_factor_integral`, `thrust_proper_velocity`, `drag_limited_proper_velocity`, `drag_limited_gamma`, `de_sitter_thrust_proper_velocity`, `de_sitter_thrust_proper_time` |
| de Sitter limit | `h_lambda`, `h_lambda_si`, `de_sitter_horizon`, `de_sitter_temperature`, `cmb_temperature`, `cmb_crossover_scale_factor`, `cmb_crossover_time`, `sink_temperature`, `lambda_tidal_acceleration`, `lambda_tidal_acceleration_si`, `recession_speed`, `de_sitter_recession_speed`, `last_light_gain`, `scale_factor_overflow_time` |

## Conventions

- `t` is cosmic time since the Big Bang of the model. "Now" is `age()`; the scale factor is 1 now.
- A cosmology carries its unit of time. Times are in that unit, lengths are light-travel times in it (`c = 1`),
  rates are per unit of time, accelerations are rapidity per unit of time. Temperatures are in kelvin, masses
  in kilograms, and functions ending in `_si` take and return SI values.
- Comoving distances `χ` are proper distances today.
- Functions of `t` return NaN for `t < 0`. Functions that can have no answer return `Option`: `None` means "never
  arrives", "was never in contact" or "outside the domain".

## The formulas

With `a_eq = (Ω_m/Ω_Λ)^(1/3)`, `H_Λ = H0 √Ω_Λ`, `x = (Ω_Λ/Ω_m) a³`, `Ω_Λ(a) = x/(1+x)` and `Ω_m(a) = 1/(1+x)`:

```text
a(t)       = a_eq sinh^(2/3)((3/2) H_Λ t)
t(a)       = 2/(3 H_Λ) · asinh((a/a_eq)^(3/2))
H(t)       = H_Λ coth((3/2) H_Λ t)                    H(a) = H0 √(Ω_m a⁻³ + Ω_Λ)

η(a)       = ∫ dt/a       = 2/(a H) · 2F1(1, 1/2; 7/6; Ω_Λ(a))          summed for x ≤ 1
η_∞ − η(a) = ∫ₜ^∞ dt/a    = 1/(a H) · 2F1(1, 1/2; 4/3; Ω_m(a))          summed for x ≥ 1
η_∞        = Γ(1/6) Γ(1/3) / (3 √π) / (a_eq H_Λ)    = 2.804 364 210 650 908 5 / (a_eq H_Λ)

∫₀ᵗ a dt   = (2/5) (a/H) · 2F1(1, 1/2; 11/6; Ω_Λ(a))                    summed for x ≤ 1
           = (a/H) · 2F1(1, 1/2; 2/3; Ω_m(a)) − K a_eq/H_Λ              summed for x ≥ 1
K          = −(2/5) Γ(11/6) Γ(−1/3) / √π             = 0.862 369 853 076 596 8
```

Each `2F1(1, 1/2; c; y)` is what the textbook form (`2F1(1/6, 1/2; 7/6; −x)` for `η`, and so on) becomes after
a Pfaff transformation, `2F1(a, b; c; z) = (1 − z)^(−b) 2F1(c − a, b; c; z/(z − 1))`. The series argument `y` is
`Ω_Λ(a)` before matter-Λ equality and `Ω_m(a)` after it, so it never exceeds 1/2, every term is positive, and
60 terms summed by Horner's rule reach the last bit. In the first form `η = 2/(aH)` is the matter-era answer and
the series is its correction; in the second `η_∞ − η = 1/(aH)` is the de Sitter answer.

Everything else is built on these:

```text
comoving distance between epochs     η(t₂) − η(t₁)
particle horizon (comoving, proper)  η,  a η
event horizon (comoving, proper)     η_∞ − η,  a (η_∞ − η)   → 1/H_Λ
reach of a signal sent at t          η_∞ − η(t)
reach mass                           c³ F³ / (2 G H_Λ),  F = a_eq H_Λ (η_∞ − η)
last contact with a target at χ      the root of η_∞ − η(t) = χ
light cone                           η_received = η_emitted + χ;  never arrives if that reaches η_∞
free particle                        γβ ∝ 1/a
rocket, proper acceleration α        d(γβ)/dt = α − H γβ   ⇒   γβ(t) = [a₀ γβ₀ + α ∫ a dt] / a(t)  → α/H_Λ
de Sitter temperature                T_dS = ħ H_Λ / (2π k_B)
microwave background                 T = T₀ / a;  T = T_dS at a = T₀ / T_dS
Λ tidal acceleration                 H_Λ² r
```

In the de Sitter limit a rocket's velocity and clock are elementary, with `U = α/H_Λ`, `γ_U = √(1 + U²)`,
`φ = asinh(γβ)`:

```text
γβ(t) = U + (γβ₀ − U) e^(−H_Λ t)
τ(t)  = 1/(H_Λ γ_U) · [ H_Λ t + 2 ln cosh((φ_U + φ(t))/2) − 2 ln cosh((φ_U + φ₀)/2) ]
      = 2/(H_Λ γ_U) · atanh( γ_U (1 − e^(−H_Λ t)) / (γ(t) + γ₀ e^(−H_Λ t)) )
```

The inverses that have no closed form (`time_at_conformal_time`, `time_at_conformal_time_remaining`, and through
them last contact and light arrival) use Newton's method on `ln η` or `ln(η_∞ − η)` as a function of `ln a`.
Both are monotone and concave, and the matter-era or de Sitter estimate starts on the side from which Newton
cannot overshoot. The worst case found takes 7 steps; the cap is 40.

One identity fell out of the tests: the particle horizon equals the event horizon when `a = a_eq/2`, that is
`η(a_eq/2) = η_∞/2`, in every flat matter + Λ universe.

## Out of scope

- **Radiation.** The closed forms are for matter and Λ. Radiation shortens the age today by 0.04% and the
  particle horizon today by 1.9%, and changes everything that looks forward from today by less than one part in
  10⁵. Do not use this crate before `a ≈ 0.01`.
- Curvature, massive neutrinos, evolving dark energy.
- Luminosity and angular-diameter distances: astropy and CCL do observational cosmology.
- The photon-count definition of Last Light (it depends on the telescope; only its aperture scaling,
  `last_light_gain`, is here). `last_light_time` is the horizon definition.
- Proper time of a rocket in full ΛCDM (no closed form; the de Sitter limit is given).

## How it was checked

`cargo test -p srkit-flrw` runs:

- **The physics oracle** (`tests/oracle.rs`): the cosmology rows of `oracle/oracle.csv` (SR-20,
  SR-21, SR-68, DT-15 to DT-18, DT-23, DT-25 to DT-29), at each row's tolerance, with the oracle's own
  `1/H_Λ` = 17.53 Gyr and with the Planck 2018 preset.
- **Golden values** (`tests/golden.rs`): 2,618 values in `tests/golden.csv`, written by `tools/crosscheck.py`
  from 50-digit mpmath quadrature of the Friedmann integrals (no hypergeometric functions on that side) and from
  astropy's `FlatLambdaCDM` with radiation off. Rows with ordinary conditioning are held to 1e-12; the worst
  measured is 4.3e-14. Rows where the question itself loses digits (two epochs 0.01% apart, a target at
  99.9999% of the horizon, `H_Λ t` in the hundreds) are held to a tolerance derived from that loss.
- **The repository's Python scripts** (`tests/scripts.rs`): the printed numbers of `lcdm_reach.py`,
  `event_horizon_windows.py` and `reach_mass_closed_form.py`. Three numbers of `lcdm_reach.py` include
  radiation and are not reproduced: age 13.79 Gyr (here 13.80), particle horizon 46.2 Gly (here 47.1) and the
  reachable share 4.6% (here 4.3%).
- **Properties** (`tests/properties.rs`): round trips, monotone and bounded conformal time, the Friedmann
  equation, the matter-era and de Sitter limits, causal light cones, last contact as the edge of the horizon,
  `γβ ∝ 1/a`, the rocket's equation of motion, additivity of proper time, the flat-space limit. 512 cases each
  from proptest's deterministic generator.
- **Bit-exact vectors** (`tests/bits.rs`): 40 results pinned to the bit. Native x86-64 (debug and release) and
  `wasm32-unknown-unknown` under Node 22 give the same bits.

To regenerate or verify the reference values (needs mpmath; astropy optional):

```sh
python tools/crosscheck.py           # rewrite tests/golden.csv
python tools/crosscheck.py --check   # fail if the checked-in file differs from a fresh run
```

## Accuracy

Against mpmath, relative error is a few units in the last place for the background and for conformal time, and
grows where the mathematics makes it grow:

- anything proportional to a power of `a` at late times carries `H_Λ t` times the rounding of `t` (9e-15 at
  `H_Λ t` = 24, 1.5e-13 at 600);
- a difference of two nearby epochs loses `η/|Δη|`;
- a signal to a target at the edge of the horizon loses `(η_∞ − η)/(η_∞ − η − χ)`;
- a burn much shorter than the Hubble time loses `1/(H Δt)` in `thrust_proper_velocity` (use
  `de_sitter_thrust_proper_velocity`, or flat-space kinematics, for those).

The scale factor overflows binary64 at `H_Λ t` = 710.5 (1.2e13 yr); `ln_scale_factor`, `hubble`,
`event_horizon_proper` and the `time_at_ln_scale_factor` inverse carry on to any finite time. The conformal time
left underflows to zero at `H_Λ t` = 745.

## References

Planck Collaboration VI 2020, A&A 641, A6 (arXiv:1807.06209); Fixsen 2009, ApJ 707, 916; Rindler 1956, MNRAS
116, 662; Davis & Lineweaver 2004 (arXiv:astro-ph/0310808); Hogg 1999 (arXiv:astro-ph/9905116); Baes, Camps &
Van De Putte 2017, MNRAS 468, 927; Loeb 2002, Phys. Rev. D 65, 047301; Busha et al. 2003, ApJ 596, 713; Gibbons
& Hawking 1977, Phys. Rev. D 15, 2738; DLMF chapter 15 (15.8.1 and 15.8.2). In the game's repository:
`docs/engine/spacetime-package-plan.md`, `docs/engine/deep-time-engine.md` section 4.
