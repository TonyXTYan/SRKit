//! Prints the bit patterns of a fixed list of results, one per line. `tools/bits.js` runs the same list through
//! the WebAssembly build and compares: the native build and the wasm32 build must agree to the last bit.
//!
//!     cargo run -q -p srkit-wasm --example bits
use srkit_wasm::*;

fn seg(values: &[f64]) {
    for i in 0..21 {
        srkit_set(i, values.get(i as usize).copied().unwrap_or(0.0));
    }
}

fn main() {
    let g = unit_constant(14);
    let mut out: Vec<f64> = vec![g];
    // a hyperbolic boost, a brake, a circle, a steered burn and a mass-burning burn
    let kinds: [&[f64]; 5] = [
        &[
            2.0,
            2000.0,
            f64::INFINITY,
            0.1,
            -0.2,
            0.3,
            0.6,
            0.0,
            0.8,
            g,
            0.0,
        ],
        &[
            2.0, 2000.0, 2010.0, 0.1, -0.2, 0.3, 0.6, 0.0, 0.8, -0.37, 4.5,
        ],
        &[
            3.0,
            2000.0,
            f64::INFINITY,
            1.0,
            2.0,
            3.0,
            7.5,
            0.11,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.7,
        ],
        &[
            5.0,
            2000.0,
            f64::INFINITY,
            0.0,
            0.0,
            0.0,
            1.2,
            -0.4,
            0.3,
            0.2,
            0.9,
            -0.1,
            0.83,
        ],
        &[
            6.0,
            2000.0,
            f64::INFINITY,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.3,
            0.71,
            0.4,
        ],
    ];
    for k in kinds {
        seg(k);
        for t in [2000.001, 2001.7, 2009.3] {
            out.push(seg_dtau(t));
            seg_at(t);
            out.extend((0..6).map(|i: u32| srkit_get(i)));
            out.push(seg_t_of_tau(seg_dtau(t)));
        }
        out.push(seg_recv(1999.0, 3.0, -2.0, 1.0, 2000.0, f64::INFINITY));
        out.push(seg_send(
            2030.0,
            3.0,
            -2.0,
            1.0,
            f64::NEG_INFINITY,
            f64::INFINITY,
        ));
    }
    // rides, reach, optics, rockets, Kepler
    for (d, cap, stops, rt) in [
        (4.2465, f64::INFINITY, 0.0, 0.0),
        (1e5, 1e4, 1.0, 1.0),
        (2.5e6, f64::INFINITY, 0.0, 0.0),
    ] {
        trip_plan(d, g, cap, stops, rt);
        out.extend((5..21).map(|i: u32| srkit_get(i)));
        let tau = srkit_get(6) * 0.37;
        trip_state(d, g, cap, stops, rt, tau, 0);
        out.extend((1..10).map(|i: u32| srkit_get(i)));
        out.push(trip_tau_of(d, g, cap, stops, rt, d * 0.5, 1));
        out.push(trip_light_cone(
            d,
            g,
            cap,
            stops,
            rt,
            d * 0.4,
            d * 0.3,
            d * 0.01,
            0,
        ));
    }
    reach_point(1.3, g, 9.0, 0.8);
    out.extend((0..2).map(|i: u32| srkit_get(i)));
    observe(3.3, 0.2, -0.5, 0.7, 1);
    out.extend((0..6).map(|i: u32| srkit_get(i)));
    boost(0.3, -1.1, 2.0, 5.0, 1.0, 2.0, 3.0);
    out.extend((0..4).map(|i: u32| srkit_get(i)));
    out.extend([
        rocket_mass_ratio(3.7, 0.3),
        sail_rapidity(0.42),
        kepler_eccentric_anomaly(2.2, 0.6),
        black_body_luminosity(6.957e8, 5772.0),
    ]);
    // cosmology and gravity
    let (h0, om) = (flrw_planck_2018(0), flrw_planck_2018(1));
    for (code, a, b) in [
        (0, 0.0, 0.0),
        (1, 5e9, 0.0),
        (4, 13.8e9, 0.0),
        (5, 0.0, 0.0),
        (6, 13.8e9, 0.0),
        (11, 1e10, 0.0),
        (13, 13.8e9, 5e9),
        (20, g, 0.0),
        (23, 0.0, 0.0),
    ] {
        out.push(flrw(code, h0, om, a, b, 0.0));
    }
    for (code, a) in [
        (0, 0.0),
        (0, 0.9),
        (0, -0.5),
        (1, 0.9),
        (2, 0.9),
        (3, 1e9),
        (4, 1e9),
        (12, 0.7),
    ] {
        out.push(gr(code, a, 0.0, 0.0, 0.0));
    }
    out.push(gr(9, 1.98841e30, 1.496e11, 2.279e11, 6.957e8));
    for v in out {
        println!("{:016x}", v.to_bits());
    }
}
