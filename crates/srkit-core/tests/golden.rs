//! Golden vectors: every case in `tests/golden/core.txt` was computed from its textbook formula at 60 digits by
//! `tools/golden_core.py`, independently of this crate. The tolerances below are the crate's accuracy promise for
//! each family, in relative terms; a looser tolerance carries its reason.
use srkit_core::{
    kepler, kinematics, optics, reach,
    segment::{Burn, Segment, Steered},
    trip::{FlipBurn, Trip},
    Vec3,
};
use std::collections::BTreeMap;

const GOLDEN: &str = include_str!("golden/core.txt");

fn parse(list: &str) -> Vec<f64> {
    list.split(',').map(|v| v.parse().unwrap()).collect()
}

/// Relative error, with `floor` as the smallest size the comparison is relative to.
fn rel(got: f64, want: f64, floor: f64) -> f64 {
    (got - want).abs() / want.abs().max(floor)
}

#[test]
fn golden_vectors() {
    let mut worst: BTreeMap<&str, f64> = BTreeMap::new();
    let mut count = 0;
    let mut failures: Vec<String> = Vec::new();
    for line in GOLDEN.lines().filter(|l| !l.is_empty()) {
        let mut parts = line.split(';');
        let case = parts.next().unwrap();
        let i = parse(parts.next().unwrap());
        let want = parse(parts.next().unwrap());
        let (got, floors, tol): (Vec<f64>, Vec<f64>, f64) = match case {
            "hyp_at" => {
                let g = Segment::hyperbolic(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, i[1], i[0]);
                let s = g.at(i[2]);
                // While braking, coordinate time resolves the remaining rapidity only to eps * sinh(phi0), so the
                // speed and proper time are compared against the starting scale.
                (
                    vec![s.x.x, s.v.x, g.dtau(i[2])],
                    vec![0.0, 1e-3, 0.0],
                    2e-13,
                )
            }
            "hyp_tau" => {
                let g = Segment::hyperbolic(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, i[1], i[0]);
                let s = g.at_tau(i[2]);
                // the product alpha * s is rounded once, and sinh amplifies that by the rapidity (up to 40 here)
                (vec![s.t, s.x.x, s.phi], vec![0.0, 0.0, 0.0], 2e-14)
            }
            "flipburn" => {
                let f = FlipBurn::new(i[0], i[1]).unwrap();
                (
                    vec![f.t, f.tau, f.phi_peak, f.gamma_peak, f.lag],
                    vec![0.0; 5],
                    4e-15,
                )
            }
            "trip" => {
                let t = Trip::new(i[0], i[1], i[2], i[3] as u32, i[4] != 0.0).unwrap();
                (
                    vec![t.tau, t.t_end, t.du, t.dv, t.gamma_peak, t.lag],
                    vec![0.0; 6],
                    // a capped leg's coast is the difference of the leg and its two burns
                    if i[2].is_finite() { 1e-11 } else { 1e-14 },
                )
            }
            "trip_state" => {
                let t = Trip::new(i[0], i[1], i[2], i[3] as u32, i[4] != 0.0).unwrap();
                let s = t.state(i[5]);
                (
                    vec![s.u, s.v, s.phi],
                    vec![0.0, 0.0, 1e-3],
                    if i[2].is_finite() { 1e-10 } else { 1e-12 },
                )
            }
            "observe" => {
                let (th, _) = (i[1], ());
                let s = optics::observe(i[0], th.cos(), th.sin(), 0.0, true).unwrap();
                (vec![s.doppler, s.theta_seen], vec![0.0, 0.0], 1e-14)
            }
            "boost" => {
                let (t, x) = kinematics::boost_rapidity(
                    Vec3::new(i[0], i[1], i[2]),
                    i[3],
                    Vec3::new(i[4], i[5], i[6]),
                );
                // compared against the size of the event, since a component can cancel to nearly nothing
                let size = (i[3].abs() + i[4].abs() + i[5].abs() + i[6].abs())
                    * Vec3::new(i[0], i[1], i[2]).norm().cosh();
                (vec![t, x.x, x.y, x.z], vec![size; 4], 4e-15)
            }
            "steered" => {
                let g = Segment::Steered(
                    Steered::new(
                        0.0,
                        f64::INFINITY,
                        Vec3::ZERO,
                        Vec3::new(i[0], i[1], i[2]),
                        Vec3::new(i[3], i[4], i[5]),
                        i[6],
                    )
                    .unwrap(),
                );
                let s = g.at(i[7]);
                let size = want[0].abs().max(want[1].abs()).max(want[2].abs());
                (
                    vec![s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z, g.dtau(i[7])],
                    vec![size, size, size, 1.0, 1.0, 1.0, 0.0],
                    // the input is a velocity: at 1 - 1e-6 of c its gamma is known to 1e-10
                    1e-9,
                )
            }
            "burn" => {
                let g =
                    Burn::new(0.0, f64::INFINITY, Vec3::ZERO, Vec3::X, i[0], i[1], i[2]).unwrap();
                let seg = Segment::Burn(g);
                let s = seg.at_tau(i[3]);
                (
                    vec![s.t, s.x.x, s.phi, g.mass_fraction(i[3])],
                    vec![0.0, s.t.abs(), 1e-3, 0.0],
                    1e-13,
                )
            }
            "kepler" => (
                vec![kepler::eccentric_anomaly(i[0], i[1])],
                vec![1.0],
                1e-14,
            ),
            "reach" => {
                let r = reach::reach(i[0], i[1], i[2], i[3]).unwrap().unwrap();
                let size = want[0].abs().max(want[1].abs());
                (r.to_vec(), vec![size, size], 1e-10)
            }
            other => panic!("unknown golden case {other}"),
        };
        assert_eq!(got.len(), want.len(), "{case}");
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            let e = rel(*g, *w, floors[k]);
            let entry = worst.entry(case).or_insert(0.0);
            if e > *entry {
                *entry = e;
            }
            if e > tol {
                failures.push(format!(
                    "{case} output {k}: got {g:e}, want {w:e}, rel {e:e} > {tol:e}\n  inputs {i:?}"
                ));
            }
        }
        count += 1;
    }
    assert_eq!(count, 610);
    assert!(
        failures.is_empty(),
        "{} golden failures:\n{}",
        failures.len(),
        failures[..failures.len().min(12)].join("\n")
    );
    for (case, e) in &worst {
        println!("{case:12} worst relative error {e:.2e}");
    }
}
