use pyo3::prelude::*;
use std::collections::HashMap;

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

impl Segment {
    /// Computes the physical length of the segment (Fast, Zero-copy, Rust-only)
    pub fn length(&self, nodes: &[Node]) -> f64 {
        let start = &nodes[self.start_idx];
        let end = &nodes[self.end_idx];
        let dx = end.x - start.x;
        let dy = end.y - start.y;
        let dz = end.z - start.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
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
    
    #[pyo3(name = "length")]
    pub fn py_length(&self, nodes: Vec<Node>) -> f64 {
        self.length(&nodes)
    }
}

/// Represents a PWS Dipole formed by two adjacent segments.
#[pyclass(from_py_object)]
#[derive(Clone, Debug)]
pub struct Dipole {
    #[pyo3(get)]
    pub seg1_idx: usize,
    #[pyo3(get)]
    pub seg2_idx: usize,
    #[pyo3(get)]
    pub junction_idx: usize,
    #[pyo3(get)]
    pub mbc_offset: f64,
}

#[pymethods]
impl Dipole {
    #[new]
    pub fn new(seg1_idx: usize, seg2_idx: usize, junction_idx: usize, mbc_offset: f64) -> Self {
        Dipole { seg1_idx, seg2_idx, junction_idx, mbc_offset }
    }
}

/// A container to hold the entire mesh and eventually compute junctions.
#[pyclass]
pub struct Mesh {
    #[pyo3(get)]
    pub nodes: Vec<Node>,
    #[pyo3(get)]
    pub segments: Vec<Segment>,
    #[pyo3(get)]
    pub dipoles: Vec<Dipole>,
}

#[pymethods]
impl Mesh {
    #[new]
    pub fn new() -> Self {
        Mesh {
            nodes: Vec::new(),
            segments: Vec::new(),
            dipoles: Vec::new(),
        }
    }

    pub fn add_node(&mut self, node: Node) -> usize {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    pub fn add_segment(&mut self, segment: Segment) {
        self.segments.push(segment);
    }

    /// Automatically scans the topology and builds N-1 dipoles for every N-way junction.
    pub fn build_dipoles(&mut self) {
        self.dipoles.clear();
        
        // Map each node index to a list of connected segment indices
        let mut connections: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, seg) in self.segments.iter().enumerate() {
            connections.entry(seg.start_idx).or_default().push(i);
            connections.entry(seg.end_idx).or_default().push(i);
        }

        // Generate the dipoles
        for (&node_idx, segs) in connections.iter() {
            if segs.len() >= 2 {
                let base_seg_idx = segs[0];
                let r_base = self.segments[base_seg_idx].radius;
                
                // Form a dipole between the base segment and every other connected segment
                for &other_seg_idx in segs.iter().skip(1) {
                    let r_other = self.segments[other_seg_idx].radius;
                    
                    // MBC Rule: Offset is the maximum of the two segment radii
                    let offset = r_base.max(r_other);
                    
                    self.dipoles.push(Dipole::new(base_seg_idx, other_seg_idx, node_idx, offset));
                }
            }
        }
    }
}

/// A Python module implemented in Rust using the modern Bound API.
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let geometry_module = PyModule::new(py, "geometry")?;

    geometry_module.add_class::<Node>()?;
    geometry_module.add_class::<Segment>()?;
    geometry_module.add_class::<Dipole>()?;
    geometry_module.add_class::<Mesh>()?;

    // m.add_function(wrap_pyfunction!(test_func, &itur_module)?)?;
    m.add_submodule(&geometry_module)?;

    let sys = py.import("sys")?;
    let sys_modules = sys.getattr("modules")?;
    sys_modules.set_item("mbc_mom.geometry", &geometry_module)?;

    Ok(())
}