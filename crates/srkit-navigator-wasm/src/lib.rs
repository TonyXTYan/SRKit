//! Small synchronous ABI: input/output arrays are fixed and owned by Rust.
//! No serialization, JavaScript callbacks or engine state cross the boundary.
// Fixed no_mangle symbols form this isolated module's linker ABI. There are
// no unsafe blocks or raw pointer dereferences; srkit-core forbids unsafe code.
use srkit_core::{
    kinematics::Vec3,
    navigator::{self, Orbit, Trip},
    Error,
};
use std::cell::RefCell;

thread_local! {
    static INPUT: RefCell<[f64; 32]> = const { RefCell::new([0.0; 32]) };
    static OUTPUT: RefCell<[f64; 32]> = const { RefCell::new([0.0; 32]) };
}

#[no_mangle]
pub extern "C" fn srkit_abi_version() -> u32 {
    1
}
#[no_mangle]
pub extern "C" fn srkit_input_ptr() -> usize {
    INPUT.with(|v| v.borrow().as_ptr() as usize)
}
#[no_mangle]
pub extern "C" fn srkit_output_ptr() -> usize {
    OUTPUT.with(|v| v.borrow().as_ptr() as usize)
}
#[no_mangle]
pub extern "C" fn srkit_set(index: u32, value: f64) {
    INPUT.with(|v| {
        if let Some(slot) = v.borrow_mut().get_mut(index as usize) {
            *slot = value;
        }
    });
}
#[no_mangle]
pub extern "C" fn srkit_get(index: u32) -> f64 {
    OUTPUT.with(|v| v.borrow().get(index as usize).copied().unwrap_or(f64::NAN))
}

fn trip(a: &[f64; 32]) -> Result<Trip, Error> {
    if !a[3].is_finite()
        || a[3] < 0.0
        || a[3].fract() != 0.0
        || a[3] > 10000.0
        || !matches!(a[4], 0.0 | 1.0)
    {
        return Err(Error::InvalidInput);
    }
    Trip::new(a[0], a[1], a[2], a[3] as u32, a[4] == 1.0)
}
fn vector(a: &[f64; 32], start: usize) -> Vec3 {
    Vec3::new(a[start], a[start + 1], a[start + 2])
}
fn execute(op: u32, a: &[f64; 32], out: &mut [f64; 32]) -> Result<(), Error> {
    match op {
        0 => {
            let p = trip(a)?;
            out[..22].copy_from_slice(&[
                p.d,
                p.a,
                p.n as f64,
                p.n_out as f64,
                if p.rt { 1.0 } else { 0.0 },
                p.leg_tau,
                p.tau,
                p.du,
                p.dv,
                p.u_end,
                p.v_end,
                p.t_end,
                p.t_out,
                p.x_end,
                p.lag,
                p.g_peak,
                p.phc,
                p.leg.a,
                p.leg.phc,
                p.leg.ta,
                p.leg.tc,
                0.0,
            ]);
        }
        1 => {
            let s = trip(a)?.state(a[5])?;
            out[..10]
                .copy_from_slice(&[s.tau, s.u, s.v, s.t, s.x, s.phi, s.dir, s.gamma, s.dh, s.dd]);
        }
        2 => {
            let p = trip(a)?;
            out[0] = if a[6] == 0.0 {
                p.tau_u(a[5])?
            } else {
                p.tau_v(a[5])?
            };
        }
        3 => {
            let p = Orbit::new(a[0], a[1])?;
            out[..8].copy_from_slice(&[
                p.a, p.gamma, p.beta, p.phi, p.lap_t, p.lap_tau, p.omega, p.sagnac,
            ]);
        }
        4 => {
            let p = Orbit::new(a[0], a[1])?;
            let s = p.state(a[2], a[3])?;
            out[..11].copy_from_slice(&[
                s.t,
                s.position.x,
                s.position.y,
                s.position.z,
                s.velocity.x,
                s.velocity.y,
                s.velocity.z,
                s.heading.x,
                s.heading.y,
                s.heading.z,
                s.angle,
            ]);
        }
        5 => {
            let v = navigator::observe(a[0], a[1], a[2], a[3], a[4] != 0.0)?;
            out[..6].copy_from_slice(&v);
        }
        6 => {
            let (t, x) = navigator::boost(vector(a, 0), a[3], vector(a, 4))?;
            out[..4].copy_from_slice(&[x.x, x.y, x.z, t]);
        }
        7 => {
            out[0] = trip(a)?.intersection(a[5], a[6], a[7], a[8] != 0.0)?;
        }
        8 => {
            out[0] = Orbit::new(a[0], a[1])?.intersection(a[2], a[3], vector(a, 4), a[7] != 0.0)?;
        }
        9 => {
            match navigator::reach(a[0], a[1], a[2], a[3])? {
                Some(v) => out[..2].copy_from_slice(&v),
                None => return Err(Error::NoIntersection),
            };
        }
        10 => {
            out[0] = navigator::traveller_emission(a[0], vector(a, 1), a[4], vector(a, 5))?;
        }
        _ => return Err(Error::InvalidInput),
    }
    if out.iter().any(|v| !v.is_finite()) {
        return Err(Error::Overflow);
    }
    Ok(())
}

/// 0 success, 1 no intersection, 2 invalid input, 3 unrepresentable output.
#[no_mangle]
pub extern "C" fn srkit_call(op: u32) -> u32 {
    let input = INPUT.with(|v| *v.borrow());
    OUTPUT.with(|v| {
        let mut out = v.borrow_mut();
        out.fill(0.0);
        match execute(op, &input, &mut out) {
            Ok(()) => 0,
            Err(Error::NoIntersection) => 1,
            Err(Error::Overflow | Error::PrecisionLoss) => 3,
            Err(_) => 2,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn abi_produces_trip_and_rejects_invalid_op() {
        for (i, v) in [4.2465, 1.0, f64::INFINITY, 0.0, 0.0]
            .into_iter()
            .enumerate()
        {
            srkit_set(i as u32, v);
        }
        assert_eq!(srkit_call(0), 0);
        assert!((srkit_get(6) - 3.542).abs() < 0.001);
        assert_eq!(srkit_call(100), 2);
        srkit_set(3, 0.5);
        assert_eq!(srkit_call(0), 2);
    }
}
