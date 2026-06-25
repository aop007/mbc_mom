use pyo3::prelude::*;

pub mod geometry;
// pub mod expj;

/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    geometry::register_module(m)?;
    // expj::register_module(m)?;
    Ok(())
}