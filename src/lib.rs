use pyo3::prelude::*;

/// Represents a 3D coordinate in space.
#[pyclass(from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub z: f64,
}

#[pymethods]
impl Node {
    #[new]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Node { x, y, z }
    }
}

/// Represents a wire segment connecting two nodes.
#[pyclass(from_py_object)]
#[derive(Clone, Debug)]
pub struct Segment {
    #[pyo3(get, set)]
    pub start_idx: usize,
    #[pyo3(get, set)]
    pub end_idx: usize,
    #[pyo3(get, set)]
    pub radius: f64,
}

#[pymethods]
impl Segment {
    #[new]
    pub fn new(start_idx: usize, end_idx: usize, radius: f64) -> Self {
        Segment {
            start_idx,
            end_idx,
            radius,
        }
    }
    
    /// Computes the physical length of the segment given the node list.
    pub fn length(&self, nodes: Vec<Node>) -> f64 {
        let start = &nodes[self.start_idx];
        let end = &nodes[self.end_idx];
        let dx = end.x - start.x;
        let dy = end.y - start.y;
        let dz = end.z - start.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// A container to hold the entire mesh and eventually compute junctions.
#[pyclass]
pub struct Mesh {
    #[pyo3(get)]
    pub nodes: Vec<Node>,
    #[pyo3(get)]
    pub segments: Vec<Segment>,
}

#[pymethods]
impl Mesh {
    #[new]
    pub fn new() -> Self {
        Mesh {
            nodes: Vec::new(),
            segments: Vec::new(),
        }
    }

    pub fn add_node(&mut self, node: Node) -> usize {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    pub fn add_segment(&mut self, segment: Segment) {
        self.segments.push(segment);
    }
}

/// A Python module implemented in Rust using the modern Bound API.
#[pymodule]
fn mbc_lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Node>()?;
    m.add_class::<Segment>()?;
    m.add_class::<Mesh>()?;
    Ok(())
}