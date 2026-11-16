//! The numbers printed by the repository's Python reference scripts in `oracle/scripts/`, reproduced
//! from the closed forms: `lcdm_reach.py`, `event_horizon_windows.py` and `reach_mass_closed_form.py`.
//!
//! The scripts integrate numerically on a grid. Since R71.10 (16 Nov 2026) they use the official constants and the
//! game's cosmology, Planck 2018 (`H0` = 67.4 km/s/Mpc, `Ω_m` = 0.315); before that they used rounded constants
//! and 67.7, 0.31. Each test compares at the digits the script prints.
//!
//! Three of `lcdm_reach.py`'s numbers cannot be reproduced by this crate, because the script includes radiation
//! (`Or = 9e-5`) and the closed forms do not: the age, the particle horizon and the share that divides by it.
//! `radiation_is_the_only_difference` pins down what the crate gives instead.

use srkit_flrw::Cosmology;

/// The scripts' constants: exact `c`, the Julian gigayear, the light-gigayear and the megaparsec.
const C: f64 = 299_792_458.0;
const GYR: f64 = 3.15576e16;
const GLY: f64 = C * GYR;
const MPC: f64 = srkit_core::units::PARSEC * 1e6;

/// `H0 = 67.4 km/s/Mpc`, `Ω_m = 0.315`, in seconds and light-seconds, as both horizon scripts set it up.
fn script_cosmology() -> Cosmology {
    Cosmology::per_second(67.4e3 / MPC, 0.315).unwrap()
}

/// Light-seconds to the script's billions of light-years.
fn gly(light_seconds: f64) -> f64 {
    light_seconds * C / GLY
}

/// The script's billions of light-years to light-seconds.
fn from_gly(g: f64) -> f64 {
    g * GLY / C
}

/// `value` printed as C's `%.{decimals}f`.
fn fixed(value: f64, decimals: usize) -> String {
    format!("{value:.decimals$}")
}

#[test]
fn lcdm_reach_py() {
    let u = script_cosmology();
    let now = u.age();
    let horizon = u.event_horizon(now);

    // "event horizon 16.7 Gly"
    assert_eq!(fixed(gly(horizon), 1), "16.7");

    // "launch at +N Gyr: reach R Gly comoving (P% of today's reach volume)"
    for (later, reach, share) in [
        (1.0, "15.71", "83.6"),
        (10.0, "9.30", "17.3"),
        (50.0, "0.95", "0.0"),
        (100.0, "0.05", "0.0"),
    ] {
        let t = now + later * GYR;
        assert_eq!(fixed(gly(u.reach(t)), 2), reach, "launch at +{later} Gyr");
        assert_eq!(
            fixed(100.0 * u.reach_volume_ratio(t, now), 1),
            share,
            "launch at +{later} Gyr"
        );
    }

    // "target X Gly comoving: arrival at +T Gyr" (printed with %.3g)
    for (target, arrival, decimals) in [
        (0.054, "0.0541", 4),
        (1.0, "1.04", 2),
        (5.0, "6.06", 2),
        (10.0, "15.8", 1),
        (15.0, "40.0", 1),
        (16.0, "55.8", 1),
    ] {
        let t = u.light_arrival_time(now, from_gly(target)).unwrap();
        assert_eq!(
            fixed((t - now) / GYR, decimals),
            arrival,
            "target {target} Gly"
        );
    }
    // 16.7 Gly is the horizon: a target at 17 Gly is never reached
    assert!(u.light_arrival_time(now, from_gly(17.0)).is_none());
}

#[test]
fn radiation_is_the_only_difference() {
    let u = script_cosmology();
    let now = u.age();
    // lcdm_reach.py prints "age 13.79 Gyr; particle horizon 46.1 Gly" and "reachable share ... 4.7%" with
    // radiation. Matter and Λ alone give 13.80 Gyr, 47.0 Gly and 4.5%.
    assert_eq!(fixed(now / GYR, 2), "13.80");
    assert_eq!(fixed(gly(u.particle_horizon(now)), 1), "47.0");
    assert_eq!(
        fixed(100.0 * u.reachable_share_of_observable(now), 1),
        "4.5"
    );
    // the script's model, with radiation (9e-5, taken out of Λ to stay flat), by 30-digit quadrature: 13.79081 Gyr,
    // 46.1427 Gly, 16.67933 Gly. The event horizon barely moves.
    let age_shift = now / GYR / 13.79081 - 1.0;
    let horizon_shift = gly(u.particle_horizon(now)) / 46.1427 - 1.0;
    let event_shift = gly(u.event_horizon(now)) / 16.67933 - 1.0;
    assert!((age_shift - 3.9e-4).abs() < 1e-5, "age shift {age_shift}");
    assert!(
        (horizon_shift - 0.0187).abs() < 1e-4,
        "particle horizon shift {horizon_shift}"
    );
    assert!(
        event_shift.abs() < 1e-4,
        "event horizon shift {event_shift}"
    );
}

#[test]
fn event_horizon_windows_py() {
    let u = script_cosmology();
    let now = u.age();
    // "event horizon now: 16.7 Gly"
    assert_eq!(fixed(gly(u.event_horizon(now)), 1), "16.7");
    // "Virgo cluster: signals sent later than 100 Gyr from now never arrive"; Coma: 69 Gyr
    for (name, megaparsecs, window) in [
        ("Virgo cluster", 16.5, "100"),
        ("Coma cluster", 100.0, "69"),
    ] {
        let chi = megaparsecs * MPC / C;
        let last = u.last_contact_time(chi).unwrap();
        assert_eq!(fixed((last - now) / GYR, 0), window, "{name}");
        // and the window is real: just inside it a signal arrives, just outside it never does
        assert!(u.light_arrival_time(last - 0.01 * GYR, chi).is_some());
        assert!(u.light_arrival_time(last + 0.01 * GYR, chi).is_none());
    }
}

#[test]
fn reach_mass_closed_form_py() {
    // the script fixes 1/H_Λ = 17.53 Gyr; the reach mass does not depend on Ω_m
    let h_lambda = 1.0 / (17.53 * GYR);
    let u = Cosmology::per_second(h_lambda / (1.0f64 - 0.315).sqrt(), 0.315).unwrap();
    // F(τ) = a_eq H_Λ (η_∞ − η): "F(0+) ~ 2.7815, closed form 2.8044"
    let f = |t: f64| u.event_horizon(t) * u.scale_factor_at_equality() * u.h_lambda();
    assert_eq!(fixed(f(0.0), 4), "2.8044");
    assert_eq!(fixed(f(1e-6 / h_lambda), 4), "2.7815");
    assert_eq!(
        u.conformal_time_limit() * u.scale_factor_at_equality() * u.h_lambda(),
        f(0.0)
    );

    // "t = T Gyr: numeric M Msun, closed form M Msun". Since R71.10 the script uses the official c, G and solar
    // mass, the same as the crate, so the crate's own value is the script's to rounding.
    let (g, m_sun) = (srkit_core::units::G_NEWTON, srkit_core::units::M_SUN);
    for (gyr, printed) in [(1.0, 3.994e23), (13.79, 2.228e22), (50.0, 4.319e19)] {
        let t = gyr * GYR;
        let script_mass = C.powi(3) / (2.0 * g * h_lambda) * f(t).powi(3) / m_sun;
        let rounded: f64 = format!("{script_mass:.3e}").parse().unwrap();
        assert_eq!(rounded, printed, "t = {gyr} Gyr");
        let crate_mass = u.reach_mass(t) / m_sun;
        assert!((crate_mass / script_mass - 1.0).abs() < 1e-13);
    }
}
