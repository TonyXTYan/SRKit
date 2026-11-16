//! # srkit-flrw
//!
//! Flat ΛCDM cosmology in closed form, as seen by an observer inside the universe.
//!
//! A [`Cosmology`] is a spatially flat universe of matter and a cosmological constant. Its scale factor is
//!
//! ```text
//! a(t) = (Ω_m/Ω_Λ)^(1/3) sinh^(2/3)((3/2) √Ω_Λ H0 t)
//! ```
//!
//! and everything else here follows from it without a numerical integration:
//!
//! * **Background.** [`Cosmology::scale_factor`], its inverse, the Hubble rate, the age and redshifts.
//! * **Conformal time.** `η = ∫ dt/a` and its finite limit `η_∞`, through Gauss hypergeometric functions
//!   summed after a Pfaff transformation ([`Cosmology::conformal_time`],
//!   [`Cosmology::conformal_time_remaining`], [`Cosmology::conformal_time_limit`]).
//! * **Horizons and reach.** The particle horizon, the cosmic event horizon, how far light sent at a given
//!   time can ever get, and the mass inside that reach ([`Cosmology::event_horizon`],
//!   [`Cosmology::reach_mass`]).
//! * **Last contact.** The latest time a signal can be sent to a comoving target and still arrive
//!   ([`Cosmology::last_contact_time`]).
//! * **Light cones.** Flat FLRW is conformally flat, so a light cone is `η_received = η_emitted + χ`
//!   ([`Cosmology::light_reception`], [`Cosmology::light_arrival_time`],
//!   [`Cosmology::light_emission_time`]).
//! * **Hubble drag.** Peculiar momentum decays as `1/a`; a rocket at constant proper acceleration saturates at
//!   `γβ = α/H_Λ` ([`Cosmology::coasting_proper_velocity`], [`Cosmology::thrust_proper_velocity`],
//!   [`Cosmology::de_sitter_thrust_proper_time`]).
//! * **The de Sitter limit.** `H_Λ`, the de Sitter horizon and its Gibbons-Hawking temperature, the cooling
//!   microwave background, Λ's tidal acceleration ([`Cosmology::de_sitter_temperature`],
//!   [`Cosmology::cmb_crossover_time`], [`Cosmology::lambda_tidal_acceleration_si`]).
//!
//! ## Conventions
//!
//! * `t` is cosmic time since the Big Bang of the model; "now" is [`Cosmology::age`]. The scale factor is 1
//!   now.
//! * A cosmology carries its unit of time, and lengths are light-travel times in that unit (`c = 1`):
//!   [`Cosmology::PLANCK_2018`] works in Julian years and light-years, [`Cosmology::PLANCK_2018_SI`] in seconds
//!   and light-seconds.
//! * Comoving distances `χ` are proper distances today.
//! * Functions of `t` return NaN for `t < 0`; functions that can have no answer return `Option`.
//!
//! ## Out of scope
//!
//! Radiation, curvature, massive neutrinos and evolving dark energy. Radiation changes the early universe (the
//! age by 0.04%, the particle horizon today by 1.9%) and leaves the future untouched to one part in 10⁵; see
//! [`Cosmology`]. Observational distance measures (luminosity and angular-diameter distances) are left to
//! astropy and CCL.
//!
//! ## Numerical contract
//!
//! `no_std`, no allocation, no `unsafe`, no threads. Every transcendental goes through the pure-Rust `libm`
//! crate, so a given build gives the same bits on every IEEE 754 binary64 platform. The hypergeometric series
//! are summed with a number of terms fixed by the argument alone, and the two root-finders
//! ([`Cosmology::time_at_conformal_time`], [`Cosmology::time_at_conformal_time_remaining`]) take a bounded
//! number of Newton steps that depends only on the arguments. Against 50-digit mpmath quadrature the relative
//! error is below 1e-13 wherever the question is well conditioned; each function says where it is not (two
//! epochs close together, a target at the very edge of the horizon, times beyond a few hundred Hubble times).
//!
//! ```
//! use srkit_flrw::Cosmology;
//!
//! let u = Cosmology::PLANCK_2018; // Julian years and light-years
//! let now = u.age();
//!
//! // The Virgo cluster, 16.5 Mpc away: 53.8 million light-years, comoving.
//! let virgo = 16.5e6 * srkit_core::units::PARSEC / srkit_core::units::LIGHT_YEAR;
//! // A signal sent now arrives in 53.9 million years ...
//! let arrival = u.light_arrival_time(now, virgo).unwrap();
//! assert!(((arrival - now) / 1e6 - 53.92).abs() < 0.01);
//! // ... and a signal sent more than 100 billion years from now never arrives.
//! let last = u.last_contact_time(virgo).unwrap();
//! assert!(((last - now) / 1e9 - 100.3).abs() < 0.1);
//! // A comoving target 20 billion light-years away is already out of reach.
//! assert!(u.light_arrival_time(now, 20e9).is_none());
//! ```
#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod conformal;
mod cosmology;
mod desitter;
mod drag;
mod series;
mod universal;

pub use conformal::Reception;
pub use cosmology::{redshift_at_scale_factor, scale_factor_at_redshift, Cosmology, KM_S_MPC};
pub use desitter::{cmb_temperature_at_scale_factor, T_CMB0};
