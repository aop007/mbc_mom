use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::geometry::{Mesh, Node, Segment};
use crate::expj;

/// Calculates the dense [Z] impedance matrix in parallel
pub fn compute_z_matrix(mesh: &Mesh, frequency_hz: f64) -> Vec<Complex64> {
    let n = mesh.dipoles.len();
    let mut z_matrix = vec![Complex64::new(0.0, 0.0); n * n];

    // Free-space wavenumber k = 2 * pi * f / c
    let k = (2.0 * PI * frequency_hz) / 299_792_458.0;

    z_matrix.par_iter_mut().enumerate().for_each(|(idx, z_val)| {
        let i = idx / n;
        let j = idx % n;

        if i <= j {
            let dip_i = &mesh.dipoles[i];
            let dip_j = &mesh.dipoles[j];

            let seg1_i = &mesh.segments[dip_i.seg1_idx];
            let seg2_i = &mesh.segments[dip_i.seg2_idx];
            let seg1_j = &mesh.segments[dip_j.seg1_idx];
            let seg2_j = &mesh.segments[dip_j.seg2_idx];

            let junc_i_seg1_is_end = seg1_i.end_idx == dip_i.junction_idx;
            let junc_i_seg2_is_end = seg2_i.end_idx == dip_i.junction_idx;
            
            let junc_j_seg1_is_end = seg1_j.end_idx == dip_j.junction_idx;
            let junc_j_seg2_is_end = seg2_j.end_idx == dip_j.junction_idx;

            // KCL Reference Directions: Current flows from seg1 -> junction -> seg2
            // seg1 should flow TOWARDS the junction
            let sign1_i = if junc_i_seg1_is_end { 1.0 } else { -1.0 };
            // seg2 should flow AWAY FROM the junction
            let sign2_i = if !junc_i_seg2_is_end { 1.0 } else { -1.0 };

            let sign1_j = if junc_j_seg1_is_end { 1.0 } else { -1.0 };
            let sign2_j = if !junc_j_seg2_is_end { 1.0 } else { -1.0 };

            let z13 = monopole_mutual(seg1_i, seg1_j, &mesh.nodes, k, dip_i.mbc_offset, dip_j.mbc_offset, junc_i_seg1_is_end, junc_j_seg1_is_end);
            let z14 = monopole_mutual(seg1_i, seg2_j, &mesh.nodes, k, dip_i.mbc_offset, dip_j.mbc_offset, junc_i_seg1_is_end, junc_j_seg2_is_end);
            let z23 = monopole_mutual(seg2_i, seg1_j, &mesh.nodes, k, dip_i.mbc_offset, dip_j.mbc_offset, junc_i_seg2_is_end, junc_j_seg1_is_end);
            let z24 = monopole_mutual(seg2_i, seg2_j, &mesh.nodes, k, dip_i.mbc_offset, dip_j.mbc_offset, junc_i_seg2_is_end, junc_j_seg2_is_end);

            if false {
                println!("[{},{}] z13: {} Ohms", i, j, z13);
                println!("[{},{}] z14: {} Ohms", i, j, z14);
                println!("[{},{}] z23: {} Ohms", i, j, z23);
                println!("[{},{}] z24: {} Ohms", i, j, z24);
            }

            // Apply the directional signs to the matrix sum
            *z_val = Complex64::new(sign1_i * sign1_j, 0.0) * z13 +
                     Complex64::new(sign1_i * sign2_j, 0.0) * z14 +
                     Complex64::new(sign2_i * sign1_j, 0.0) * z23 +
                     Complex64::new(sign2_i * sign2_j, 0.0) * z24;
        }
    });

    // Mirror the upper triangle to the lower triangle (Exact Reciprocity)
    for i in 0..n {
        for j in 0..i {
            z_matrix[i * n + j] = z_matrix[j * n + i];
        }
    }

    z_matrix
}

/// Robust 2D Numerical Quadrature for Monopole Mutual Impedance.
/// Safely bypasses the 1/R singularity due to the MBC multiradius offset.
fn monopole_mutual(
    seg_a: &Segment, seg_b: &Segment,
    nodes: &[Node], k: f64,
    offset_a: f64, offset_b: f64,
    junc_a_is_end: bool, junc_b_is_end: bool
) -> Complex64 {
    let n_a1 = &nodes[seg_a.start_idx];
    let n_a2 = &nodes[seg_a.end_idx];
    let n_b1 = &nodes[seg_b.start_idx];
    let n_b2 = &nodes[seg_b.end_idx];

    let len_a = seg_a.length(nodes);
    let len_b = seg_b.length(nodes);

    let s_hat = [
        (n_a2.x - n_a1.x) / len_a,
        (n_a2.y - n_a1.y) / len_a,
        (n_a2.z - n_a1.z) / len_a,
    ];

    let t_hat = [
        (n_b2.x - n_b1.x) / len_b,
        (n_b2.y - n_b1.y) / len_b,
        (n_b2.z - n_b1.z) / len_b,
    ];

    let s_dot_t = s_hat[0]*t_hat[0] + s_hat[1]*t_hat[1] + s_hat[2]*t_hat[2];

    // The MBC multiradius rule bounds the Green's function
    let mbc_offset = offset_a.max(offset_b);
    let offset_sq = mbc_offset * mbc_offset;

    // --- NEW: Dynamic Quadrature Resolution ---
    // A thin wire (large aspect ratio) needs extremely fine resolution 
    // to capture the sharp 1/R near-field singularity peak.
    let aspect_ratio = len_a.max(len_b) / mbc_offset;
    
    // Scale steps dynamically. Minimum 50 for thick wires, capped at 2000 
    // for ultra-thin wires to maintain sub-second performance in Rust.
    let steps = (aspect_ratio * 0.8) as usize;
    let steps = steps.clamp(50, 2000);

    let ds = len_a / (steps as f64);
    let dt = len_b / (steps as f64);

    let sin_kl_a = (k * len_a).sin();
    let sin_kl_b = (k * len_b).sin();

    let mut sum = Complex64::new(0.0, 0.0);
    let j_cplx = Complex64::new(0.0, 1.0);

    for i in 0..steps {
        let s = (i as f64 + 0.5) * ds;
        
        // PWS Current Shape & Derivative for Segment A
        let (i_a, di_a) = if junc_a_is_end {
            ( (k * s).sin() / sin_kl_a, k * (k * s).cos() / sin_kl_a )
        } else {
            ( (k * (len_a - s)).sin() / sin_kl_a, -k * (k * (len_a - s)).cos() / sin_kl_a )
        };

        let px = n_a1.x + s_hat[0] * s;
        let py = n_a1.y + s_hat[1] * s;
        let pz = n_a1.z + s_hat[2] * s;

        for j in 0..steps {
            let t = (j as f64 + 0.5) * dt;

            // PWS Current Shape & Derivative for Segment B
            let (i_b, di_b) = if junc_b_is_end {
                ( (k * t).sin() / sin_kl_b, k * (k * t).cos() / sin_kl_b )
            } else {
                ( (k * (len_b - t)).sin() / sin_kl_b, -k * (k * (len_b - t)).cos() / sin_kl_b )
            };

            let qx = n_b1.x + t_hat[0] * t;
            let qy = n_b1.y + t_hat[1] * t;
            let qz = n_b1.z + t_hat[2] * t;

            let dx = px - qx;
            let dy = py - qy;
            let dz = pz - qz;
            
            // 3D Distance with MBC bounded offset
            let dist = (dx*dx + dy*dy + dz*dz + offset_sq).sqrt();

            let green = (-j_cplx * k * dist).exp() / dist;

            // Symmetric Mutual Impedance Integrand (Tilston & Balmain 1988, Eq 5)
            let term = s_dot_t * i_a * i_b - (di_a * di_b) / (k * k);
            
            sum += green * term;
        }
    }

    // Multiply by (j * \eta * k) / (4 * \pi) * ds * dt
    let multiplier = j_cplx * (376.7303 / (4.0 * PI)) * k;
    sum * multiplier * ds * dt
}
