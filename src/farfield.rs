use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::geometry::{Mesh, Segment};

/// Computes the Radiation Intensity U(theta, phi) over a grid of angles.
pub fn compute_pattern(
    mesh: &Mesh,
    currents: Vec<Complex64>,
    freq_hz: f64,
    thetas: Vec<f64>,
    phis: Vec<f64>,
) -> Vec<f64> {
    let k = (2.0 * PI * freq_hz) / 299_792_458.0;
    let eta = 376.7303; // Free-space impedance

    // 1. Precompute the exact vector moment of every radiating segment
    struct Radiator {
        center: [f64; 3],
        moment: [Complex64; 3],
    }
    let mut radiators = Vec::new();

    for (i, dip) in mesh.dipoles.iter().enumerate() {
        let i_dip = currents[i];
        let seg1 = &mesh.segments[dip.seg1_idx];
        let seg2 = &mesh.segments[dip.seg2_idx];

        let junc_seg1_is_end = seg1.end_idx == dip.junction_idx;
        let junc_seg2_is_end = seg2.end_idx == dip.junction_idx;

        // KCL Flow exactly as defined in impedance.rs
        let sign1 = if junc_seg1_is_end { 1.0 } else { -1.0 };
        let sign2 = if !junc_seg2_is_end { 1.0 } else { -1.0 };

        let mut process_seg = |seg: &Segment, sign: f64| {
            let start = &mesh.nodes[seg.start_idx];
            let end = &mesh.nodes[seg.end_idx];
            
            // Midpoint phase center approximation for electrically short segments
            let center = [
                (start.x + end.x) / 2.0,
                (start.y + end.y) / 2.0,
                (start.z + end.z) / 2.0,
            ];
            
            let len = seg.length(&mesh.nodes);
            let u = [
                (end.x - start.x) / len,
                (end.y - start.y) / len,
                (end.z - start.z) / len,
            ];

            // Exact integral of the PWS current shape over the segment
            let kl = k * len;
            let q = if kl < 1e-5 {
                len / 2.0 // Taylor expansion limit for extreme thin/short wires
            } else {
                (1.0 - kl.cos()) / (k * kl.sin())
            };

            let moment_mag = i_dip * Complex64::new(sign * q, 0.0);
            let moment = [
                moment_mag * u[0],
                moment_mag * u[1],
                moment_mag * u[2],
            ];

            radiators.push(Radiator { center, moment });
        };

        process_seg(seg1, sign1);
        process_seg(seg2, sign2);
    }

    let mut u_grid = vec![0.0; thetas.len() * phis.len()];

    // 2. Compute Far-Field Grid in Parallel
    u_grid.par_iter_mut().enumerate().for_each(|(idx, u_val)| {
        let t_idx = idx / phis.len();
        let p_idx = idx % phis.len();
        let theta = thetas[t_idx];
        let phi = phis[p_idx];

        let st = theta.sin();
        let ct = theta.cos();
        let sp = phi.sin();
        let cp = phi.cos();

        // Spherical unit vectors mapped to Cartesian
        let r_hat = [st * cp, st * sp, ct];
        let theta_hat = [ct * cp, ct * sp, -st];
        let phi_hat = [-sp, cp, 0.0];

        let mut n_vec = [Complex64::new(0.0, 0.0); 3];

        // Sum the phase-shifted moments
        for rad in &radiators {
            let phase = k * (r_hat[0] * rad.center[0] + r_hat[1] * rad.center[1] + r_hat[2] * rad.center[2]);
            let phase_cplx = Complex64::new(0.0, phase).exp();

            n_vec[0] += rad.moment[0] * phase_cplx;
            n_vec[1] += rad.moment[1] * phase_cplx;
            n_vec[2] += rad.moment[2] * phase_cplx;
        }

        // Project onto transverse observation plane
        let n_theta = n_vec[0] * theta_hat[0] + n_vec[1] * theta_hat[1] + n_vec[2] * theta_hat[2];
        let n_phi = n_vec[0] * phi_hat[0] + n_vec[1] * phi_hat[1] + n_vec[2] * phi_hat[2];

        let n_trans_sq = n_theta.norm_sqr() + n_phi.norm_sqr();
        
        // Radiation Intensity U(theta, phi) in Watts/steradian
        *u_val = (k * k * eta / (32.0 * PI * PI)) * n_trans_sq;
    });

    u_grid
}