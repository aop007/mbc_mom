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
    #[pyo3(get)]
    pub is_monopole: bool,
}

#[pymethods]
impl Dipole {
    #[new]
    pub fn new(seg1_idx: usize, seg2_idx: usize, junction_idx: usize, mbc_offset: f64, is_monopole: bool) -> Self {
        Dipole { seg1_idx, seg2_idx, junction_idx, mbc_offset, is_monopole }
    }
}

#[pyclass(from_py_object)]
#[derive(Clone, Debug)]
pub struct GroundPlane {
    #[pyo3(get)]
    pub is_pec: bool,
    #[pyo3(get)]
    pub sigma: f64,
    #[pyo3(get)]
    pub eps_r: f64,
    #[pyo3(get)]
    pub use_sommerfeld: bool,
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
    #[pyo3(get)]
    pub ground_plane: Option<GroundPlane>,
}

#[pymethods]
impl Mesh {
    #[new]
    pub fn new() -> Self {
        Mesh {
            nodes: Vec::new(),
            segments: Vec::new(),
            dipoles: Vec::new(),
            ground_plane: None,
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
        let mut dipoles = Vec::new();
        let has_ground = self.ground_plane.is_some();
        
        // 1. Standard Two-Segment Dipoles
        for i in 0..self.segments.len() {
            for j in (i + 1)..self.segments.len() {
                let s1 = &self.segments[i];
                let s2 = &self.segments[j];

                let junction = if s1.start_idx == s2.start_idx || s1.start_idx == s2.end_idx {
                    Some(s1.start_idx)
                } else if s1.end_idx == s2.start_idx || s1.end_idx == s2.end_idx {
                    Some(s1.end_idx)
                } else {
                    None
                };

                if let Some(j_idx) = junction {
                    let mbc_offset = s1.radius.max(s2.radius);
                    dipoles.push(Dipole {
                        seg1_idx: i,
                        seg2_idx: j,
                        junction_idx: j_idx,
                        mbc_offset,
                        is_monopole: false,
                    });
                }
            }
        }

        // 2. Monopole Dipoles (Grounded Segments)
        if has_ground {
            for i in 0..self.segments.len() {
                let s = &self.segments[i];
                let n_start = &self.nodes[s.start_idx];
                let n_end = &self.nodes[s.end_idx];

                // If a node touches the z=0 plane (within a tiny floating-point tolerance)
                let tol = 1e-12;
                let (is_grounded, j_idx) = if n_start.z.abs() < tol {
                    (true, s.start_idx)
                } else if n_end.z.abs() < tol {
                    (true, s.end_idx)
                } else {
                    (false, 0)
                };

                if is_grounded {
                    dipoles.push(Dipole {
                        seg1_idx: i,
                        seg2_idx: i, // Placeholder: Monopoles only have one physical segment
                        junction_idx: j_idx,
                        mbc_offset: s.radius, // Image has the exact same radius
                        is_monopole: true,
                    });
                }
            }
        }

        self.dipoles = dipoles;
    }

    /// Sets the environment to a Perfectly Electric Conducting (PEC) half-space
    pub fn set_pec_ground(&mut self) {
        self.ground_plane = Some(GroundPlane {
            is_pec: true,
            sigma: 0.0,
            eps_r: 1.0,
            use_sommerfeld: false,
        });
    }

    /// Sets the environment to a Homogeneous Dielectric half-space
    pub fn set_real_ground(&mut self, sigma: f64, eps_r: f64, use_sommerfeld: bool) {
        self.ground_plane = Some(GroundPlane {
            is_pec: false,
            sigma,
            eps_r,
            use_sommerfeld,
        });
    }
}

/// A Python module implemented in Rust using the modern Bound API.
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let geometry_module = PyModule::new(py, "geometry")?;

    geometry_module.add_class::<Node>()?;
    geometry_module.add_class::<Segment>()?;
    geometry_module.add_class::<Dipole>()?;
    geometry_module.add_class::<GroundPlane>()?;
    geometry_module.add_class::<Mesh>()?;

    // m.add_function(wrap_pyfunction!(test_func, &itur_module)?)?;
    m.add_submodule(&geometry_module)?;

    let sys = py.import("sys")?;
    let sys_modules = sys.getattr("modules")?;
    sys_modules.set_item("mbc_mom.geometry", &geometry_module)?;

    Ok(())
}