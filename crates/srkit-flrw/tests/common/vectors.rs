//! The bit-exact test vectors: one call of every function family, at an early epoch, today and a late epoch.
//!
//! This file uses nothing but `srkit_flrw` and `core`, so the same list can be compiled for another target
//! (it has been run as WebAssembly under Node) and compared bit for bit with `tests/bits.rs`.

use srkit_flrw::{Cosmology, T_CMB0};

/// Number of vectors.
pub const COUNT: usize = 40;

/// Name and value of every vector, for the Planck 2018 preset in Julian years and light-years.
pub fn vectors() -> [(&'static str, f64); COUNT] {
    let u = Cosmology::PLANCK_2018;
    let now = u.age();
    let early = 3.0e5;
    let late = 2.0e11;
    let far = 9.0e12;
    let g = u.acceleration_from_si(9.806_65);
    let chi = 5.0e9;
    let some = |value: Option<f64>| value.unwrap_or(f64::NAN);
    [
        ("age", now),
        ("h_lambda", u.h_lambda()),
        ("scale_factor_at_equality", u.scale_factor_at_equality()),
        ("scale_factor early", u.scale_factor(early)),
        ("scale_factor late", u.scale_factor(late)),
        ("ln_scale_factor far", u.ln_scale_factor(far)),
        ("time_at_scale_factor", u.time_at_scale_factor(0.25)),
        ("time_at_ln_scale_factor", u.time_at_ln_scale_factor(900.0)),
        ("hubble early", u.hubble(early)),
        ("hubble_at_scale_factor", u.hubble_at_scale_factor(0.5)),
        ("lookback_time", u.lookback_time(1.0)),
        ("redshift_between", u.redshift_between(early, now)),
        ("conformal_time early", u.conformal_time(early)),
        ("conformal_time now", u.conformal_time(now)),
        ("conformal_time_limit", u.conformal_time_limit()),
        (
            "conformal_time_remaining now",
            u.conformal_time_remaining(now),
        ),
        (
            "conformal_time_remaining far",
            u.conformal_time_remaining(far),
        ),
        (
            "time_at_conformal_time",
            some(u.time_at_conformal_time(2.0e10)),
        ),
        (
            "time_at_conformal_time_remaining",
            some(u.time_at_conformal_time_remaining(3.0e9)),
        ),
        ("comoving_distance", u.comoving_distance(early, late)),
        (
            "comoving_distance_at_redshift",
            u.comoving_distance_at_redshift(3.0),
        ),
        ("particle_horizon_proper", u.particle_horizon_proper(now)),
        ("event_horizon_proper", u.event_horizon_proper(now)),
        (
            "reachable_share_of_observable",
            u.reachable_share_of_observable(now),
        ),
        ("reach_mass", u.reach_mass(now)),
        ("last_contact_time", some(u.last_contact_time(chi))),
        ("last_reply_time", some(u.last_reply_time(chi))),
        ("light_arrival_time", some(u.light_arrival_time(now, chi))),
        ("light_emission_time", some(u.light_emission_time(now, chi))),
        (
            "light_reception",
            some(u.light_reception(4.0e10, chi).map(|r| r.time)),
        ),
        (
            "coasting_proper_velocity",
            u.coasting_proper_velocity(3.0, now, late),
        ),
        ("coasting_rapidity", u.coasting_rapidity(25.0, now, late)),
        ("coasting_gamma", u.coasting_gamma(7.0, early, now)),
        (
            "scale_factor_integral early",
            u.scale_factor_integral(early),
        ),
        ("scale_factor_integral late", u.scale_factor_integral(late)),
        (
            "thrust_proper_velocity",
            u.thrust_proper_velocity(0.0, g, now, late),
        ),
        (
            "de_sitter_thrust_proper_time short",
            u.de_sitter_thrust_proper_time(0.0, g, 10.0),
        ),
        (
            "de_sitter_thrust_proper_time long",
            u.de_sitter_thrust_proper_time(0.0, g, 1.0e100),
        ),
        ("de_sitter_temperature", u.de_sitter_temperature()),
        ("cmb_crossover_time", u.cmb_crossover_time(T_CMB0)),
    ]
}
