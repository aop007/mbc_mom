use pyo3::prelude::*;

mod geometry;
// pub mod expj;

#[pyfunction]
pub fn test_interface() -> () {
    println!("Rust Interface Running!");
}

/// Represents a 3D coordinate in space.
#[pyclass(from_py_object)]
#[derive(Clone, Debug)]
pub struct Foo {
    #[pyo3(get, set)]
    pub bar: f64,
}

#[pymethods]
impl Foo {
    #[new]
    pub fn new(bar: f64) -> Self {
        Foo { bar }
    }
}

/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_mom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // expj::register_module(m)?;
    m.add_function(wrap_pyfunction!(test_interface, m)?)?;
    m.add_class::<Foo>()?;

    geometry::register_module(m)?;
    
    Ok(())
}