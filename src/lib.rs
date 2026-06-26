use pyo3::prelude::*;
use num_complex::Complex64;

mod geometry;
mod impedance;
mod farfield;
mod sommerfeld;
pub mod expj;

use geometry::{Mesh};

#[pyfunction]
pub fn test_interface() -> () {
    println!("Rust Interface Running!");
}

/// Python wrapper to compute the Z matrix. Returns a flat 1D list of complex numbers.
#[pyfunction]
pub fn compute_impedance_matrix(mesh: &Mesh, frequency_hz: f64) -> Vec<Complex64> {
    impedance::compute_z_matrix(mesh, frequency_hz)
}

#[pyfunction]
pub fn compute_far_field(mesh: &Mesh, currents: Vec<Complex64>, freq_hz: f64, thetas: Vec<f64>, phis: Vec<f64>) -> Vec<f64> {
    farfield::compute_pattern(mesh, currents, freq_hz, thetas, phis)
}

/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_mom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // expj::register_module(m)?;
    m.add_function(wrap_pyfunction!(test_interface, m)?)?;
    m.add_function(wrap_pyfunction!(compute_impedance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(compute_far_field, m)?)?;

    geometry::register_module(m)?;
    
    Ok(())
}