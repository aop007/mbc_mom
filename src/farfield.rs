use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::constants::{C, EPS_0, ETA};
use crate::geometry::{Mesh, Segment};

pub fn compute_pattern(
    mesh: &Mesh,
    currents: Vec<Complex64>,
    freq_hz: f64,
    thetas_phis: Vec<(f64, f64)>,
) -> Vec<f64> {
    let omega = 2.0 * PI * freq_hz;
    let k = omega / C;

    // 1. Extract Ground Parameters
    let has_ground = mesh.ground_plane.is_some();
    let (is_pec, eps_r, sigma) = if let Some(g) = &mesh.ground_plane {
        (g.is_pec, g.eps_r, g.sigma)
    } else {
        (false, 1.0, 0.0)
    };

    // Calculate Complex Permittivity
    let eps_c = Complex64::new(eps_r, -sigma / (omega * EPS_0));

    struct Radiator {
        center: [f64; 3],
        moment: [Complex64; 3],
    }
    let mut radiators = Vec::new();

    // 2. Precompute ONLY Physical Radiators
    for (i, dip) in mesh.dipoles.iter().enumerate() {
        let i_dip = currents[i];
        let seg1 = &mesh.segments[dip.seg1_idx];
        let sign1 = if seg1.end_idx == dip.junction_idx { 1.0 } else { -1.0 };

        let mut process_seg = |seg: &Segment, sign: f64| {
            let start = &mesh.nodes[seg.start_idx];
            let end = &mesh.nodes[seg.end_idx];
            
            let center = [(start.x + end.x) / 2.0, (start.y + end.y) / 2.0, (start.z + end.z) / 2.0];
            let len = seg.length(&mesh.nodes);
            let u = [(end.x - start.x) / len, (end.y - start.y) / len, (end.z - start.z) / len];

            let kl = k * len;
            let q = if kl < 1e-5 { len / 2.0 } else { (1.0 - kl.cos()) / (k * kl.sin()) };
            let moment_mag = i_dip * Complex64::new(sign * q, 0.0);
            
            radiators.push(Radiator {
                center,
                moment: [moment_mag * u[0], moment_mag * u[1], moment_mag * u[2]],
            });
        };

        process_seg(seg1, sign1);
        if !dip.is_monopole {
            let seg2 = &mesh.segments[dip.seg2_idx];
            let sign2 = if seg2.end_idx != dip.junction_idx { 1.0 } else { -1.0 };
            process_seg(seg2, sign2);
        }
    }

    let mut u_grid = vec![0.0; thetas_phis.len()];

    // 3. Compute Grid with Fresnel Reflection
    u_grid.par_iter_mut().enumerate().for_each(|(idx, u_val)| {
        let theta = thetas_phis[idx].0;
        let phi = thetas_phis[idx].1;

        // Mask underground radiation
        if has_ground && theta >= (PI / 2.0) {
            *u_val = 0.0;
            return;
        }

        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();

        let r_hat = [st * cp, st * sp, ct];
        let theta_hat = [ct * cp, ct * sp, -st];
        let phi_hat = [-sp, cp, 0.0];

        // Evaluate Fresnel Coefficients for this specific elevation angle
        let (gamma_v, gamma_h) = if has_ground {
            if is_pec {
                (Complex64::new(1.0, 0.0), Complex64::new(-1.0, 0.0))
            } else {
                let cost = Complex64::new(ct, 0.0);
                let sint_sq = Complex64::new(st * st, 0.0);
                let root = (eps_c - sint_sq).sqrt();
                
                let g_v = (eps_c * cost - root) / (eps_c * cost + root);
                let g_h = (cost - root) / (cost + root);
                (g_v, g_h)
            }
        } else {
            (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0))
        };

        let mut n_theta_total = Complex64::new(0.0, 0.0);
        let mut n_phi_total = Complex64::new(0.0, 0.0);

        for rad in &radiators {
            // Direct Ray
            let phase_dir = k * (r_hat[0] * rad.center[0] + r_hat[1] * rad.center[1] + r_hat[2] * rad.center[2]);
            let exp_dir = Complex64::new(0.0, phase_dir).exp();
            
            let e_theta_dir = (rad.moment[0] * theta_hat[0] + rad.moment[1] * theta_hat[1] + rad.moment[2] * theta_hat[2]) * exp_dir;
            let e_phi_dir = (rad.moment[0] * phi_hat[0] + rad.moment[1] * phi_hat[1] + rad.moment[2] * phi_hat[2]) * exp_dir;

            // Reflected Ray (Phase shifted to the underground image location)
            let mut e_theta_ref = Complex64::new(0.0, 0.0);
            let mut e_phi_ref = Complex64::new(0.0, 0.0);

            if has_ground {
                let phase_ref = k * (r_hat[0] * rad.center[0] + r_hat[1] * rad.center[1] + r_hat[2] * (-rad.center[2]));
                let exp_ref = Complex64::new(0.0, phase_ref).exp();
                
                e_theta_ref = (rad.moment[0] * theta_hat[0] + rad.moment[1] * theta_hat[1] + rad.moment[2] * theta_hat[2]) * exp_ref;
                e_phi_ref = (rad.moment[0] * phi_hat[0] + rad.moment[1] * phi_hat[1] + rad.moment[2] * phi_hat[2]) * exp_ref;
            }

            // Superposition: Direct + (Reflected * Gamma)
            n_theta_total += e_theta_dir + (e_theta_ref * gamma_v);
            n_phi_total += e_phi_dir + (e_phi_ref * gamma_h);
        }

        *u_val = (k * k * ETA / (32.0 * PI * PI)) * (n_theta_total.norm_sqr() + n_phi_total.norm_sqr());
    });

    u_grid
}