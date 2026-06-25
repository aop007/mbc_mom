use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

use crate::geometry::{Mesh, Node, Segment};

/// Calculates the dense [Z] impedance matrix in parallel
pub fn compute_z_matrix(mesh: &Mesh, frequency_hz: f64) -> Vec<Complex64> {
    let n = mesh.dipoles.len();
    let mut z_matrix = vec![Complex64::new(0.0, 0.0); n * n];

    let omega = 2.0 * PI * frequency_hz;
    let mu_0 = 4.0 * PI * 1e-7;
    let eps_0 = 8.854187817e-12;
    let k = Complex64::new(omega * (mu_0 * eps_0).sqrt(), 0.0); // Free-space wavenumber

    // Parallelize across the flat N x N matrix
    z_matrix.par_iter_mut().enumerate().for_each(|(idx, z_val)| {
        let i = idx / n;
        let j = idx % n;

        // MBC enforces exact reciprocity, so Z_ij == Z_ji. 
        // We only compute the upper triangle and diagonal to save 50% compute time!
        if i <= j {
            let dipole_i = &mesh.dipoles[i];
            let dipole_j = &mesh.dipoles[j];

            // A dipole is composed of two segments (monopoles).
            // Let Dipole I = Monopoles 1 & 2. Dipole J = Monopoles 3 & 4.
            // Z_ij = Z_13 + Z_14 + Z_23 + Z_24
            
            let z13 = monopole_mutual(
                &mesh.segments[dipole_i.seg1_idx], 
                &mesh.segments[dipole_j.seg1_idx], 
                &mesh.nodes, k, dipole_i.mbc_offset, dipole_j.mbc_offset
            );
            
            let z14 = monopole_mutual(
                &mesh.segments[dipole_i.seg1_idx], 
                &mesh.segments[dipole_j.seg2_idx], 
                &mesh.nodes, k, dipole_i.mbc_offset, dipole_j.mbc_offset
            );
            
            let z23 = monopole_mutual(
                &mesh.segments[dipole_i.seg2_idx], 
                &mesh.segments[dipole_j.seg1_idx], 
                &mesh.nodes, k, dipole_i.mbc_offset, dipole_j.mbc_offset
            );
            
            let z24 = monopole_mutual(
                &mesh.segments[dipole_i.seg2_idx], 
                &mesh.segments[dipole_j.seg2_idx], 
                &mesh.nodes, k, dipole_i.mbc_offset, dipole_j.mbc_offset
            );

            *z_val = z13 + z14 + z23 + z24;
        }
    });

    // Mirror the upper triangle to the lower triangle
    for i in 0..n {
        for j in 0..i {
            z_matrix[i * n + j] = z_matrix[j * n + i];
        }
    }

    z_matrix
}

/// Computes the mutual impedance between two filamentary monopoles.
/// In the full implementation, this calls your EXPJ analytic formulations.
fn monopole_mutual(
    seg_a: &Segment, 
    seg_b: &Segment, 
    nodes: &Vec<Node>, 
    k: Complex64,
    offset_a: f64,
    offset_b: f64
) -> Complex64 {
    // ------------------------------------------------------------------
    // TODO: Phase 3.5 - Insert the Richmond/Tilston rigorous integrals here.
    // For now, we return a dummy placeholder so the architecture compiles.
    // ------------------------------------------------------------------
    Complex64::new(1.0, 0.0) 
}