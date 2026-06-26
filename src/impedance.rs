use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::geometry::{Mesh, Dipole};
use crate::expj;

struct SubSegment {
    p1: [f64; 3],
    p2: [f64; 3],
    radius: f64,
    sign: f64,
    junc_is_p2: bool,
    weight: Complex64,
}

/// Generates the mathematical sub-segments (physical and virtual) for a given dipole.
fn get_subsegments(dip: &Dipole, mesh: &Mesh, is_source: bool, gamma: Complex64) -> Vec<SubSegment> {
    let mut subs = Vec::new();
    
    // Segment 1 (Always physical)
    let seg1 = &mesh.segments[dip.seg1_idx];
    let n1_start = &mesh.nodes[seg1.start_idx];
    let n1_end = &mesh.nodes[seg1.end_idx];
    let junc_is_end1 = seg1.end_idx == dip.junction_idx;
    let sign1 = if junc_is_end1 { 1.0 } else { -1.0 };
    
    // Testing Segment: Only physical space
    subs.push(SubSegment {
        p1: [n1_start.x, n1_start.y, n1_start.z],
        p2: [n1_end.x, n1_end.y, n1_end.z],
        radius: seg1.radius,
        sign: sign1,
        junc_is_p2: junc_is_end1,
        weight: Complex64::new(1.0, 0.0),
    });

    // Source Segment: Image radiates into physical space
    if is_source && gamma != Complex64::new(0.0, 0.0) {
        subs.push(SubSegment {
            p1: [n1_start.x, n1_start.y, -n1_start.z],
            p2: [n1_end.x, n1_end.y, -n1_end.z],
            radius: seg1.radius,
            sign: sign1,
            junc_is_p2: junc_is_end1,
            weight: gamma,
        });
    }
    
    // If it's a standard dipole, process the second half
    if !dip.is_monopole {
        let seg2 = &mesh.segments[dip.seg2_idx];
        let n2_start = &mesh.nodes[seg2.start_idx];
        let n2_end = &mesh.nodes[seg2.end_idx];
        let junc_is_end2 = seg2.end_idx == dip.junction_idx;
        let sign2 = if !junc_is_end2 { 1.0 } else { -1.0 };
        
        subs.push(SubSegment {
            p1: [n2_start.x, n2_start.y, n2_start.z],
            p2: [n2_end.x, n2_end.y, n2_end.z],
            radius: seg2.radius,
            sign: sign2,
            junc_is_p2: junc_is_end2,
            weight: Complex64::new(1.0, 0.0),
        });
        
        if is_source && gamma != Complex64::new(0.0, 0.0) {
            subs.push(SubSegment {
                p1: [n2_start.x, n2_start.y, -n2_start.z],
                p2: [n2_end.x, n2_end.y, -n2_end.z],
                radius: seg2.radius,
                sign: sign2,
                junc_is_p2: junc_is_end2,
                weight: gamma,
            });
        }
    }
    
    subs
}

/// Calculates the dense [Z] impedance matrix in parallel
pub fn compute_z_matrix(mesh: &Mesh, frequency_hz: f64) -> Vec<Complex64> {
    let n = mesh.dipoles.len();
    let mut z_matrix = vec![Complex64::new(0.0, 0.0); n * n];

    // Free-space wavenumber k = 2 * pi * f / c
    let k = (2.0 * PI * frequency_hz) / 299_792_458.0;

    let has_ground = mesh.ground_plane.is_some();
    let is_pec = mesh.ground_plane.as_ref().map_or(false, |g| g.is_pec);

    let gamma = if has_ground {
        if is_pec {
            Complex64::new(-1.0, 0.0) // Perfect Electric Conductor
        } else {
            // TODO: RCA / Sommerfeld integration for Dielectric Ground.
            // Placeholder: Default to PEC for stability until RCA is implemented.
            Complex64::new(-1.0, 0.0) 
        }
    } else {
        Complex64::new(0.0, 0.0) // Free Space
    };

    z_matrix.par_iter_mut().enumerate().for_each(|(idx, z_val)| {
        let i = idx / n;
        let j = idx % n;

        if i <= j {
            let dip_i = &mesh.dipoles[i];
            let dip_j = &mesh.dipoles[j];

            // Testing segments (Physical wire only)
            let subs_i = get_subsegments(dip_i, mesh, false, gamma);
            // Source segments (Physical wire + Ground Images)
            let subs_j = get_subsegments(dip_j, mesh, true, gamma);

            let mut z_ij = Complex64::new(0.0, 0.0);

            // Dynamically calculate the cross-coupling between all physical and virtual segments
            for sub_i in &subs_i {
                for sub_j in &subs_j {
                    let mut_z = monopole_mutual(sub_i, sub_j, k);
                    
                    // Multiply by KCL signs and the image reflection weight
                    let total_weight = Complex64::new(sub_i.sign * sub_j.sign, 0.0) * sub_j.weight;
                    z_ij += total_weight * mut_z;
                }
            }

            *z_val = z_ij;
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

/// Robust 2D Numerical Quadrature mapped between two arbitrary 3D line segments
fn monopole_mutual(sub_a: &SubSegment, sub_b: &SubSegment, k: f64) -> Complex64 {
    let dx_a = sub_a.p2[0] - sub_a.p1[0];
    let dy_a = sub_a.p2[1] - sub_a.p1[1];
    let dz_a = sub_a.p2[2] - sub_a.p1[2];
    let len_a = (dx_a*dx_a + dy_a*dy_a + dz_a*dz_a).sqrt();
    let s_hat = [dx_a/len_a, dy_a/len_a, dz_a/len_a];

    let dx_b = sub_b.p2[0] - sub_b.p1[0];
    let dy_b = sub_b.p2[1] - sub_b.p1[1];
    let dz_b = sub_b.p2[2] - sub_b.p1[2];
    let len_b = (dx_b*dx_b + dy_b*dy_b + dz_b*dz_b).sqrt();
    let t_hat = [dx_b/len_b, dy_b/len_b, dz_b/len_b];

    let s_dot_t = s_hat[0]*t_hat[0] + s_hat[1]*t_hat[1] + s_hat[2]*t_hat[2];

    let mbc_offset = sub_a.radius.max(sub_b.radius);
    let offset_sq = mbc_offset * mbc_offset;

    let aspect_ratio = len_a.max(len_b) / mbc_offset;
    let steps = ((aspect_ratio * 0.8) as usize).clamp(10, 2000); 

    let ds = len_a / (steps as f64);
    let dt = len_b / (steps as f64);

    let sin_kl_a = (k * len_a).sin();
    let sin_kl_b = (k * len_b).sin();

    let mut sum = Complex64::new(0.0, 0.0);
    let j_cplx = Complex64::new(0.0, 1.0);

    for i in 0..steps {
        let s = (i as f64 + 0.5) * ds;
        let (i_a, di_a) = if sub_a.junc_is_p2 {
            ( (k * s).sin() / sin_kl_a, k * (k * s).cos() / sin_kl_a )
        } else {
            ( (k * (len_a - s)).sin() / sin_kl_a, -k * (k * (len_a - s)).cos() / sin_kl_a )
        };

        let px = sub_a.p1[0] + s_hat[0] * s;
        let py = sub_a.p1[1] + s_hat[1] * s;
        let pz = sub_a.p1[2] + s_hat[2] * s;

        for j in 0..steps {
            let t = (j as f64 + 0.5) * dt;
            let (i_b, di_b) = if sub_b.junc_is_p2 {
                ( (k * t).sin() / sin_kl_b, k * (k * t).cos() / sin_kl_b )
            } else {
                ( (k * (len_b - t)).sin() / sin_kl_b, -k * (k * (len_b - t)).cos() / sin_kl_b )
            };

            let qx = sub_b.p1[0] + t_hat[0] * t;
            let qy = sub_b.p1[1] + t_hat[1] * t;
            let qz = sub_b.p1[2] + t_hat[2] * t;

            let dist = ((px-qx).powi(2) + (py-qy).powi(2) + (pz-qz).powi(2) + offset_sq).sqrt();
            let green = (-j_cplx * k * dist).exp() / dist;
            
            sum += green * (s_dot_t * i_a * i_b - (di_a * di_b) / (k * k));
        }
    }

    sum * j_cplx * (376.7303 / (4.0 * PI)) * k * ds * dt
}
