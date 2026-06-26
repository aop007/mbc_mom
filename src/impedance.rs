use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::geometry::{Mesh, Dipole};
use crate::sommerfeld::{GroundPhysics, SommerfeldTable};

struct SubSegment {
    p1: [f64; 3],
    p2: [f64; 3],
    radius: f64,
    sign: f64,
    junc_is_p2: bool,
    is_image: bool
}

/// Generates the mathematical sub-segments (physical and virtual) for a given dipole.
fn get_subsegments(dip: &Dipole, mesh: &Mesh, is_source: bool, has_ground: bool) -> Vec<SubSegment> {
    let mut subs = Vec::new();
    
    let seg1 = &mesh.segments[dip.seg1_idx];
    let n1_start = &mesh.nodes[seg1.start_idx];
    let n1_end = &mesh.nodes[seg1.end_idx];
    let junc_is_end1 = seg1.end_idx == dip.junction_idx;
    let sign1 = if junc_is_end1 { 1.0 } else { -1.0 };
    
    subs.push(SubSegment {
        p1: [n1_start.x, n1_start.y, n1_start.z],
        p2: [n1_end.x, n1_end.y, n1_end.z],
        radius: seg1.radius,
        sign: sign1,
        junc_is_p2: junc_is_end1,
        is_image: false,
    });

    if is_source && has_ground {
        subs.push(SubSegment {
            p1: [n1_start.x, n1_start.y, -n1_start.z],
            p2: [n1_end.x, n1_end.y, -n1_end.z],
            radius: seg1.radius,
            sign: sign1,
            junc_is_p2: junc_is_end1,
            is_image: true,
        });
    }
    
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
            is_image: false,
        });
        
        if is_source && has_ground {
            subs.push(SubSegment {
                p1: [n2_start.x, n2_start.y, -n2_start.z],
                p2: [n2_end.x, n2_end.y, -n2_end.z],
                radius: seg2.radius,
                sign: sign2,
                junc_is_p2: junc_is_end2,
                is_image: true,
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
    let (is_pec, eps_r, sigma, use_sommerfeld) = if let Some(g) = &mesh.ground_plane {
        (g.is_pec, g.eps_r, g.sigma, g.use_sommerfeld)
    } else {
        (false, 1.0, 0.0, false)
    };

    let eps_0 = 8.8541878128e-12;
    let omega = 2.0 * PI * frequency_hz;
    let eps_c = Complex64::new(eps_r, -sigma / (omega * eps_0));

    // Generate the Sommerfeld LUT if requested
    let lut = if has_ground && !is_pec && use_sommerfeld {
        let physics = GroundPhysics::new(frequency_hz, eps_r, sigma);
        Some(SommerfeldTable::generate(&physics))
    } else {
        None
    };

    z_matrix.par_iter_mut().enumerate().for_each(|(idx, z_val)| {
        let i = idx / n;
        let j = idx % n;

        if i <= j {
            let dip_i = &mesh.dipoles[i];
            let dip_j = &mesh.dipoles[j];

            let subs_i = get_subsegments(dip_i, mesh, false, has_ground);
            let subs_j = get_subsegments(dip_j, mesh, true, has_ground);

            let mut z_ij = Complex64::new(0.0, 0.0);

            for sub_i in &subs_i {
                for sub_j in &subs_j {
                    let mut_z = monopole_mutual(sub_i, sub_j, k);
                    let mut total_weight = Complex64::new(sub_i.sign * sub_j.sign, 0.0);

                    // --- NEW: Near-Field RCA ---
                    if sub_j.is_image {
                        let cx_i = (sub_i.p1[0] + sub_i.p2[0]) / 2.0;
                        let cy_i = (sub_i.p1[1] + sub_i.p2[1]) / 2.0;
                        let cz_i = (sub_i.p1[2] + sub_i.p2[2]) / 2.0;

                        let cx_j = (sub_j.p1[0] + sub_j.p2[0]) / 2.0;
                        let cy_j = (sub_j.p1[1] + sub_j.p2[1]) / 2.0;
                        let cz_j = (sub_j.p1[2] + sub_j.p2[2]) / 2.0;

                        let dx = cx_i - cx_j;
                        let dy = cy_i - cy_j;
                        let dz = cz_i - cz_j; // Since j is an image, cz_j is negative
                        let r = (dx*dx + dy*dy + dz*dz).sqrt();

                        let gamma = if is_pec {
                            Complex64::new(-1.0, 0.0)
                        } else if let Some(lut_ref) = &lut {
                            // --- EXACT SOMMERFELD COUPLING ---
                            let rho = (dx*dx + dy*dy).sqrt();
                            let z_sum = cz_i.abs() + cz_j.abs();
                            
                            // 1. Get the exact numerical integral from the 2D grid
                            let exact_g = lut_ref.interpolate(rho, z_sum);
                            
                            // 2. Evaluate the ideal spatial image Green's function
                            let ideal_g = (-Complex64::new(0.0, 1.0) * k * r).exp() / r;
                            
                            // 3. Extract the Effective Reflection Coefficient
                            if ideal_g.norm() > 1e-12 {
                                -(exact_g / ideal_g) 
                            } else {
                                Complex64::new(0.0, 0.0)
                            }
                        } else {
                            // --- FALLBACK RCA COUPLING ---
                            let cost = if r == 0.0 { 1.0 } else { (dz / r).abs() };
                            let sint_sq = 1.0 - cost * cost;
                            
                            let cost_cplx = Complex64::new(cost, 0.0);
                            let sint_sq_cplx = Complex64::new(sint_sq, 0.0);
                            let root = (eps_c - sint_sq_cplx).sqrt();
                            
                            // Calculate TM (Vertical) Reflection Coefficient
                            let g_v = (eps_c * cost_cplx - root) / (eps_c * cost_cplx + root);
                            
                            // We use -g_v as the scalar to correctly map to our Z-flipped geometry.
                            // This ensures vertical currents are correctly scaled by +g_v.
                            -g_v 
                        };
                        total_weight *= gamma;
                    }

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
