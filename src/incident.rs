use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::constants::C;
use crate::geometry::Mesh;

/// Generates the flat 1D coordinate arrays to be evaluated by the Python PropagationModel.
pub fn get_incident_eval_points(
    mesh: &Mesh, 
    points_per_seg: usize
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let num_points = mesh.segments.len() * points_per_seg;
    let mut xs = Vec::with_capacity(num_points);
    let mut ys = Vec::with_capacity(num_points);
    let mut zs = Vec::with_capacity(num_points);

    for seg in &mesh.segments {
        let start = &mesh.nodes[seg.start_idx];
        let end = &mesh.nodes[seg.end_idx];
        
        let dx = end.x - start.x;
        let dy = end.y - start.y;
        let dz = end.z - start.z;
        
        let ds = 1.0 / (points_per_seg as f64);
        
        // Extract midpoint coordinates for numerical quadrature
        for i in 0..points_per_seg {
            let s_frac = (i as f64 + 0.5) * ds; 
            xs.push(start.x + dx * s_frac);
            ys.push(start.y + dy * s_frac);
            zs.push(start.z + dz * s_frac);
        }
    }

    (xs, ys, zs)
}

/// Integrates the Python-evaluated E-fields back into the MoM [V] matrix.
pub fn compute_incident_v_matrix(
    mesh: &Mesh,
    freq_hz: f64,
    ex: Vec<Complex64>,
    ey: Vec<Complex64>,
    ez: Vec<Complex64>,
    points_per_seg: usize,
) -> Vec<Complex64> {
    let k = (2.0 * PI * freq_hz) / C;
    let mut v_matrix = vec![Complex64::new(0.0, 0.0); mesh.dipoles.len()];

    v_matrix.par_iter_mut().enumerate().for_each(|(i, v_val)| {
        let dip = &mesh.dipoles[i];
        let mut v_m = Complex64::new(0.0, 0.0);

        let process_seg = |seg_idx: usize, is_seg1: bool| -> Complex64 {
            let seg = &mesh.segments[seg_idx];
            let start = &mesh.nodes[seg.start_idx];
            let end = &mesh.nodes[seg.end_idx];
            
            let dx = end.x - start.x;
            let dy = end.y - start.y;
            let dz = end.z - start.z;
            let len = (dx * dx + dy * dy + dz * dz).sqrt();
            let u_hat = if len > 1e-14 { [dx/len, dy/len, dz/len] } else { [0.0, 0.0, 1.0] };
            
            let junc_is_end = seg.end_idx == dip.junction_idx;
            let sign = if is_seg1 {
                if junc_is_end { 1.0 } else { -1.0 }
            } else {
                if !junc_is_end { 1.0 } else { -1.0 }
            };

            let ds_len = len / (points_per_seg as f64);
            let sin_kl = (k * len).sin();
            let mut seg_v = Complex64::new(0.0, 0.0);

            for step in 0..points_per_seg {
                let offset = seg_idx * points_per_seg + step;
                let e_vec = [ex[offset], ey[offset], ez[offset]];
                
                let s = (step as f64 + 0.5) * ds_len; 
                
                // Piecewise-Sinusoidal (PWS) Basis Function Weight
                let current_mag = if junc_is_end {
                    (k * s).sin() / sin_kl
                } else {
                    (k * (len - s)).sin() / sin_kl
                };
                
                // Project the incident E-field vector onto the wire's unit vector
                let e_dot_u = e_vec[0] * u_hat[0] + e_vec[1] * u_hat[1] + e_vec[2] * u_hat[2];
                
                // Superposition: (E dot U) * I(s) * ds
                seg_v += e_dot_u * Complex64::new(sign * current_mag, 0.0) * ds_len;
            }
            seg_v
        };

        v_m += process_seg(dip.seg1_idx, true);
        if !dip.is_monopole {
            v_m += process_seg(dip.seg2_idx, false);
        }

        *v_val = v_m;
    });

    v_matrix
}