//! Bit-exact vectors: the determinism contract of `docs/engine/spacetime-package-plan.md` section 3.
//!
//! Every transcendental goes through the pure-Rust `libm`, so these bits must be the same on every IEEE 754
//! binary64 platform. They were recorded on x86-64 Linux and reproduced, bit for bit, by the same code built
//! for `wasm32-unknown-unknown` and run under Node 22. A change to any of them changes what a stored game
//! replays to: treat it as a breaking change and say so in the release notes.
//!
//! To print the current bits after an intended change: `FLRW_PRINT_BITS=1 cargo test -p srkit-flrw --test bits
//! -- --nocapture`.

#[path = "common/vectors.rs"]
mod vectors;

use vectors::{vectors, COUNT};

const EXPECTED: [u64; COUNT] = [
    0x4209_b28e_aca0_0eb3, // age = 1.3796234644007177e10
    0x3dcf_5d20_db7b_075e, // h_lambda = 5.7050388228730436e-11
    0x3fe8_b319_2e3d_a34e, // scale_factor_at_equality = 7.718625929966068e-1
    0x3f46_02c8_d0cb_0eaa, // scale_factor early = 6.717186493791678e-4
    0x40e5_6bfe_e406_b02f, // scale_factor late = 4.387196533522044e4
    0x4080_05dc_0d47_be11, // ln_scale_factor far = 5.127324472050551e2
    0x41df_eb1b_6006_a0cf, // time_at_scale_factor = 2.142006656103565e9
    0x42ac_b7f0_2691_878b, // time_at_ln_scale_factor = 1.5788166826179771e13
    0x3ec2_a42f_9631_10b0, // hubble early = 2.2222222227104345e-6
    0x3de0_f5df_8cad_1faa, // hubble_at_scale_factor = 1.2340345961943753e-10
    0x41fd_9e4b_8dbf_6b71, // lookback_time = 7.950612699963731e9
    0x4097_3edf_c638_5788, // redshift_between = 1.4877185295871186e3
    0x41d3_f71c_052d_99c0, // conformal_time early = 1.3398466767125092e9
    0x4225_e396_d293_889d, // conformal_time now = 4.700631277776682e10
    0x422d_a7d1_fca5_5dac, // conformal_time_limit = 6.368480417868295e10
    0x420f_10ec_a847_543b, // conformal_time_remaining now = 1.6678491400916128e10
    0x13d3_de58_2eca_dbeb, // conformal_time_remaining far = 3.688662171849893e-213
    0x41cd_bf4c_0dca_72cf, // time_at_conformal_time = 9.981522195816287e8
    0x4224_4b56_ecdd_681b, // time_at_conformal_time_remaining = 4.358166692670333e10
    0x422d_080c_eb1e_e7e1, // comoving_distance = 6.234455796745289e10
    0x4213_c33a_29d7_40ee, // comoving_distance_at_redshift = 2.1219936885813408e10
    0x4225_e396_d293_889c, // particle_horizon_proper = 4.7006312777766815e10
    0x420f_10ec_a847_543e, // event_horizon_proper = 1.6678491400916134e10
    0x3fa6_decb_6b89_9320, // reachable_share_of_observable = 4.466853797638337e-2
    0x4add_8e36_7339_6f2a, // reach_mass = 4.4232128972544645e52
    0x4220_20ac_e04c_81f2, // last_contact_time = 3.46338386302538e10
    0x4214_fd9a_ad47_7235, // last_reply_time = 2.253852961786153e10
    0x4212_7dec_61e0_233f, // light_arrival_time = 1.985551372003442e10
    0x4201_c296_4ab8_762a, // light_emission_time = 9.534884183057697e9
    0x4206_33c3_3d7a_0c6d, // light_reception = 1.1919714223256067e10
    0x3f11_ecf5_2a20_1ef0, // coasting_proper_velocity = 6.838079801252026e-5
    0x402c_9f37_5b22_e5a9, // coasting_rapidity = 1.4310969207773114e1
    0x3ff0_000b_5adb_3e5d, // coasting_gamma = 1.0000108288840217e0
    0x405e_3a32_e736_abc8, // scale_factor_integral early = 1.2090935688342063e2
    0x4305_db27_7724_cc63, // scale_factor_integral late = 7.689921325039484e14
    0x4210_d9f1_b5dd_81fa, // thrust_proper_velocity = 1.809399333537693e10
    0x4007_7b03_ea55_a034, // de_sitter_thrust_proper_time short = 2.9350660617237967e0
    0x5291_5d0a_d86d_8f87, // de_sitter_thrust_proper_time long = 5.526557137262459e89
    0x39c6_498b_9f9d_3848, // de_sitter_temperature = 2.197696148628334e-30
    0x4271_dbc7_f58d_e773, // cmb_crossover_time = 1.2272281417904656e12
];

#[test]
fn every_vector_has_the_recorded_bits() {
    let got = vectors();
    if std::env::var_os("FLRW_PRINT_BITS").is_some() {
        for (name, value) in got {
            println!("    0x{:016x}, // {name} = {value:e}", value.to_bits());
        }
    }
    let mut wrong = Vec::new();
    for ((name, value), bits) in got.into_iter().zip(EXPECTED) {
        assert!(value.is_finite(), "{name} is not finite");
        if value.to_bits() != bits {
            wrong.push(format!(
                "{name}: 0x{:016x} ({value:e}), recorded 0x{bits:016x} ({:e})",
                value.to_bits(),
                f64::from_bits(bits)
            ));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} vectors changed:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}
