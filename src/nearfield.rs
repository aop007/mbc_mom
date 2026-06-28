use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::constants::{C, EPS_0, ETA, J};
use crate::geometry::{Mesh};
use crate::sommerfeld::{GroundPhysics, SommerfeldTable};

// A flattened physical or virtual radiator for the near-field integration
struct RadiatingSegment {
    p1: [f64; 3],
    p2: [f64; 3],
    radius: f64,
    sign: f64,
    junc_is_p2: bool,
    is_image: bool,
    i_dip: Complex64,
}

pub fn compute_grid(
    mesh: &Mesh,
    currents: Vec<Complex64>,
    freq_hz: f64,
    xs: Vec<f64>,
    ys: Vec<f64>,
    zs: Vec<f64>,
) -> (Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>) {
    
    let omega = 2.0 * PI * freq_hz;
    let k = omega / C;

    let has_ground = mesh.ground_plane.is_some();
    let (is_pec, eps_r, sigma, use_sommerfeld) = if let Some(g) = &mesh.ground_plane {
        (g.is_pec, g.eps_r, g.sigma, g.use_sommerfeld)
    } else {
        (false, 1.0, 0.0, false)
    };

    let eps_c = Complex64::new(eps_r, -sigma / (omega * EPS_0));

    // Optional: Precompute Sommerfeld Ground Physics
    let lut = if has_ground && !is_pec && use_sommerfeld {
        let physics = GroundPhysics::new(freq_hz, eps_r, sigma);
        Some(SommerfeldTable::generate(&physics))
    } else {
        None
    };

    // 1. Flatten the entire mesh into independent radiating line segments
    let mut radiators = Vec::new();
    for (i, dip) in mesh.dipoles.iter().enumerate() {
        let i_dip = currents[i];
        
        let mut extract_subs = |seg_idx: usize, is_seg1: bool| {
            let seg = &mesh.segments[seg_idx];
            let start = &mesh.nodes[seg.start_idx];
            let end = &mesh.nodes[seg.end_idx];
            let junc_is_end = seg.end_idx == dip.junction_idx;
            
            // KCL Direction
            let sign = if is_seg1 {
                if junc_is_end { 1.0 } else { -1.0 }
            } else {
                if !junc_is_end { 1.0 } else { -1.0 }
            };

            radiators.push(RadiatingSegment {
                p1: [start.x, start.y, start.z], p2: [end.x, end.y, end.z],
                radius: seg.radius, sign, junc_is_p2: junc_is_end,
                is_image: false, i_dip,
            });

            // If ground is active, build the virtual image segment
            if has_ground {
                radiators.push(RadiatingSegment {
                    p1: [start.x, start.y, -start.z], p2: [end.x, end.y, -end.z],
                    radius: seg.radius, sign, junc_is_p2: junc_is_end,
                    is_image: true, i_dip,
                });
            }
        };

        extract_subs(dip.seg1_idx, true);
        if !dip.is_monopole { extract_subs(dip.seg2_idx, false); }
    }

    // 2. Parallel Evaluation over the user's flat spatial coordinate arrays
    let points = xs.len();
    let results: Vec<_> = (0..points).into_par_iter().map(|idx| {
        let obs = [xs[idx], ys[idx], zs[idx]];

        // Do not calculate fields underneath the ground plane
        if has_ground && obs[2] < 0.0 {
            return (
                Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)
            );
        }

        let mut e_total = [Complex64::new(0.0, 0.0); 3];
        let mut h_total = [Complex64::new(0.0, 0.0); 3];

        for rad in &radiators {
            let dx_s = rad.p2[0] - rad.p1[0];
            let dy_s = rad.p2[1] - rad.p1[1];
            let dz_s = rad.p2[2] - rad.p1[2];
            let len = (dx_s*dx_s + dy_s*dy_s + dz_s*dz_s).sqrt();
            let u_hat = [dx_s/len, dy_s/len, dz_s/len];

            // Determine Ground Reflection Coefficient
            let gamma = if rad.is_image {
                let img_cz = (rad.p1[2] + rad.p2[2]) / 2.0;
                let dz = obs[2] - img_cz;
                let dx = obs[0] - (rad.p1[0] + rad.p2[0]) / 2.0;
                let dy = obs[1] - (rad.p1[1] + rad.p2[1]) / 2.0;
                let r_dist = (dx*dx + dy*dy + dz*dz).sqrt();

                if is_pec {
                    Complex64::new(-1.0, 0.0)
                } else if let Some(lut_ref) = &lut {
                    let rho = (dx*dx + dy*dy).sqrt();
                    let z_sum = obs[2].abs() + img_cz.abs();
                    let exact_g = lut_ref.interpolate(rho, z_sum);
                    let ideal_g = (-J * k * r_dist).exp() / r_dist;
                    if ideal_g.norm() > 1e-12 { -(exact_g / ideal_g) } else { Complex64::new(0.0, 0.0) }
                } else {
                    let cost = if r_dist == 0.0 { 1.0 } else { (dz / r_dist).abs() };
                    let root = (eps_c - Complex64::new(1.0 - cost*cost, 0.0)).sqrt();
                    let g_v = (eps_c * cost - root) / (eps_c * cost + root);
                    -g_v
                }
            } else {
                Complex64::new(1.0, 0.0)
            };

            // 1D Numerical Quadrature for dE and dH
            let steps = 40;
            let ds = len / (steps as f64);
            let sin_kl = (k * len).sin();

            for i in 0..steps {
                let s = (i as f64 + 0.5) * ds;
                
                // Current distribution I(s)
                let current_mag = if rad.junc_is_p2 {
                    (k * s).sin() / sin_kl
                } else {
                    (k * (len - s)).sin() / sin_kl
                };
                let current = rad.i_dip * Complex64::new(rad.sign * current_mag, 0.0) * gamma;

                // Source position along the segment
                let px = rad.p1[0] + u_hat[0] * s;
                let py = rad.p1[1] + u_hat[1] * s;
                let pz = rad.p1[2] + u_hat[2] * s;

                let dx = obs[0] - px;
                let dy = obs[1] - py;
                let dz = obs[2] - pz;
                let r_geom = (dx*dx + dy*dy + dz*dz).sqrt();
                
                // MBC Boundary: Clamp evaluation distance to the wire's surface to prevent singularity
                let r = r_geom.max(rad.radius); 
                let r_hat = if r_geom > 1e-14 { [dx/r_geom, dy/r_geom, dz/r_geom] } else { [0.0, 0.0, 1.0] };

                let exp_term = (-J * k * r).exp();
                let kr = k * r;
                let kr_sq = kr * kr;
                let u_dot_r = u_hat[0]*r_hat[0] + u_hat[1]*r_hat[1] + u_hat[2]*r_hat[2];

                // Electric Field (dE)
                let e_const = (-J * ETA / (4.0 * PI * k)) * (current * ds * exp_term / (r * r * r));
                let t1 = Complex64::new(kr_sq - 1.0, -kr);
                let t2 = Complex64::new(3.0 - kr_sq, 3.0 * kr);
                
                e_total[0] += e_const * (t1 * u_hat[0] + t2 * u_dot_r * r_hat[0]);
                e_total[1] += e_const * (t1 * u_hat[1] + t2 * u_dot_r * r_hat[1]);
                e_total[2] += e_const * (t1 * u_hat[2] + t2 * u_dot_r * r_hat[2]);

                // Magnetic Field (dH)
                let h_const = (current * ds * exp_term) / (4.0 * PI * r * r);
                let h_term = Complex64::new(1.0, kr);
                
                h_total[0] += h_const * h_term * (u_hat[1]*r_hat[2] - u_hat[2]*r_hat[1]);
                h_total[1] += h_const * h_term * (u_hat[2]*r_hat[0] - u_hat[0]*r_hat[2]);
                h_total[2] += h_const * h_term * (u_hat[0]*r_hat[1] - u_hat[1]*r_hat[0]);
            }
        }

        (e_total[0], e_total[1], e_total[2], h_total[0], h_total[1], h_total[2])
    }).collect();

    // 3. Unzip the results back into fast flat arrays for Python
    let mut ex = Vec::with_capacity(points); let mut ey = Vec::with_capacity(points); let mut ez = Vec::with_capacity(points);
    let mut hx = Vec::with_capacity(points); let mut hy = Vec::with_capacity(points); let mut hz = Vec::with_capacity(points);
    
    for res in results {
        ex.push(res.0); ey.push(res.1); ez.push(res.2);
        hx.push(res.3); hy.push(res.4); hz.push(res.5);
    }

    (ex, ey, ez, hx, hy, hz)
}