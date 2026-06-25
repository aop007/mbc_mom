use pyo3::prelude::*;

mod geometry;
pub mod expj;

#[pyfunction]
pub fn test_interface() -> () {
    println!("Rust Interface Running!");
}


/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_mom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // expj::register_module(m)?;
    m.add_function(wrap_pyfunction!(test_interface, m)?)?;

    geometry::register_module(m)?;
    
    Ok(())
}