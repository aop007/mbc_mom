use pyo3::prelude::*;
use num_complex::Complex64;

mod geometry;
mod impedance;
mod farfield;
mod sommerfeld;
mod nearfield;
mod incident;
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
pub fn compute_far_field(mesh: &Mesh, currents: Vec<Complex64>, freq_hz: f64, thetas_phis: Vec<(f64, f64)>) -> Vec<f64> {
    farfield::compute_pattern(mesh, currents, freq_hz, thetas_phis)
}

#[pyfunction]
pub fn compute_near_field(
    mesh: &Mesh, 
    currents: Vec<Complex64>, 
    freq_hz: f64, 
    xs: Vec<f64>, 
    ys: Vec<f64>, 
    zs: Vec<f64>
) -> (Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>, Vec<Complex64>) {
    nearfield::compute_grid(mesh, currents, freq_hz, xs, ys, zs)
}

#[pyfunction]
#[pyo3(signature = (mesh, points_per_seg=7))]
pub fn get_incident_eval_points(
    mesh: &Mesh, 
    points_per_seg: usize
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    incident::get_incident_eval_points(mesh, points_per_seg)
}

#[pyfunction]
#[pyo3(signature = (mesh, freq_hz, ex, ey, ez, points_per_seg=7))]
pub fn compute_incident_v_matrix(
    mesh: &Mesh,
    freq_hz: f64,
    ex: Vec<Complex64>,
    ey: Vec<Complex64>,
    ez: Vec<Complex64>,
    points_per_seg: usize,
) -> Vec<Complex64> {
    incident::compute_incident_v_matrix(mesh, freq_hz, ex, ey, ez, points_per_seg)
}

/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_mom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // expj::register_module(m)?;
    m.add_function(wrap_pyfunction!(test_interface, m)?)?;
    m.add_function(wrap_pyfunction!(compute_impedance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(compute_far_field, m)?)?;
    m.add_function(wrap_pyfunction!(compute_near_field, m)?)?;
    m.add_function(wrap_pyfunction!(get_incident_eval_points, m)?)?;
    m.add_function(wrap_pyfunction!(compute_incident_v_matrix, m)?)?;

    geometry::register_module(m)?;
    
    Ok(())
}