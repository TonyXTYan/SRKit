//! The numbers printed by the repository's Python reference scripts in `tools/physics-checks/`, reproduced
//! from the closed forms: `lcdm_reach.py`, `event_horizon_windows.py` and `reach_mass_closed_form.py`.
//!
//! The scripts integrate numerically on a grid and use rounded constants (`c` = 2.998e8 m/s, 1 Gyr = 3.156e16 s,
//! 1 Gly = 9.461e24 m, 1 Mpc = 3.0857e22 m). Each test converts with the script's own constants and compares
//! at the digits the script prints.
//!
//! Three of `lcdm_reach.py`'s numbers cannot be reproduced by this crate, because the script includes radiation
//! (`Or = 9e-5`) and the closed forms do not: the age, the particle horizon and the share that divides by it.
//! `radiation_is_the_only_difference` pins down what the crate gives instead.

use srkit_flrw::Cosmology;

/// The scripts' rounded constants.
const C: f64 = 2.998e8;
const GYR: f64 = 3.156e16;
const GLY: f64 = 9.461e24;
const MPC: f64 = 3.0857e22;

/// `H0 = 67.7 km/s/Mpc`, `Ω_m = 0.31`, in seconds and light-seconds, as both horizon scripts set it up.
fn script_cosmology() -> Cosmology {
    Cosmology::per_second(67.7 / 3.0857e19, 0.31).unwrap()
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

    // "event horizon 16.6 Gly"
    assert_eq!(fixed(gly(horizon), 1), "16.6");

    // "launch at +N Gyr: reach R Gly comoving (P% of today's reach volume)"
    for (later, reach, share) in [
        (1.0, "15.59", "83.5"),
        (10.0, "9.19", "17.1"),
        (50.0, "0.92", "0.0"),
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
        (5.0, "6.07", 2),
        (10.0, "15.8", 1),
        (15.0, "40.8", 1),
        (16.0, "58.6", 1),
    ] {
        let t = u.light_arrival_time(now, from_gly(target)).unwrap();
        assert_eq!(
            fixed((t - now) / GYR, decimals),
            arrival,
            "target {target} Gly"
        );
    }
    // 16.6 Gly is the horizon: a target at 17 Gly is never reached
    assert!(u.light_arrival_time(now, from_gly(17.0)).is_none());
}

#[test]
fn radiation_is_the_only_difference() {
    let u = script_cosmology();
    let now = u.age();
    // lcdm_reach.py prints "age 13.79 Gyr; particle horizon 46.2 Gly" and "reachable share ... 4.6%" with
    // radiation. Matter and Λ alone give 13.80 Gyr, 47.1 Gly and 4.3%.
    assert_eq!(fixed(now / GYR, 2), "13.80");
    assert_eq!(fixed(gly(u.particle_horizon(now)), 1), "47.1");
    assert_eq!(
        fixed(100.0 * u.reachable_share_of_observable(now), 1),
        "4.3"
    );
    // the script's values, with radiation, from 50-digit quadrature of the same model: 13.7903 Gyr,
    // 46.2298 Gly, 16.5608 Gly (the event horizon barely moves: 16.5609 Gly without radiation)
    let age_shift = now / GYR / 13.7903 - 1.0;
    let horizon_shift = gly(u.particle_horizon(now)) / 46.2298 - 1.0;
    let event_shift = gly(u.event_horizon(now)) / 16.5608 - 1.0;
    assert!((age_shift - 4.2e-4).abs() < 1e-5, "age shift {age_shift}");
    assert!(
        (horizon_shift - 0.0189).abs() < 1e-4,
        "particle horizon shift {horizon_shift}"
    );
    assert!(
        event_shift.abs() < 1e-5,
        "event horizon shift {event_shift}"
    );
}

#[test]
fn event_horizon_windows_py() {
    let u = script_cosmology();
    let now = u.age();
    // "event horizon now: 16.6 Gly"
    assert_eq!(fixed(gly(u.event_horizon(now)), 1), "16.6");
    // "Virgo cluster: signals sent later than 99 Gyr from now never arrive"; Coma: 68 Gyr
    for (name, megaparsecs, window) in
        [("Virgo cluster", 16.5, "99"), ("Coma cluster", 100.0, "68")]
    {
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

    // "t = T Gyr: numeric M Msun, closed form M Msun", with the script's G and solar mass
    let (g, m_sun) = (6.674e-11, 1.989e30);
    for (gyr, printed) in [(1.0, 3.994e23), (13.79, 2.227e22), (50.0, 4.318e19)] {
        let t = gyr * GYR;
        let script_mass = C.powi(3) / (2.0 * g * h_lambda) * f(t).powi(3) / m_sun;
        let rounded: f64 = format!("{script_mass:.3e}").parse().unwrap();
        assert_eq!(rounded, printed, "t = {gyr} Gyr");
        // the crate's own value uses the exact c and CODATA G in place of the script's rounded ones, which
        // lowers it by 1.2e-4: the same closed form, different constants
        let constants = (srkit_core::units::C / C).powi(3) * g / srkit_core::units::G_NEWTON;
        assert!((constants - 1.0 + 1.204e-4).abs() < 1e-7);
        let crate_mass = u.reach_mass(t) / m_sun;
        assert!((crate_mass / (script_mass * constants) - 1.0).abs() < 1e-13);
    }
}
