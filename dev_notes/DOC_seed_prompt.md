
---

### 📋 Copy-Paste This Into a New Thread

**Prompt:**

> I am developing a high-performance Method of Moments (MoM) computational electromagnetics solver. The heavy numerical physics (impedance matrix assembly, exact Sommerfeld integrals, parallelized near/far field evaluation) are written in Rust using `rayon`. This Rust backend is compiled into a Python package called `mbc_mom` using PyO3 and Maturin.
> I need to write the documentation for this package. To give you the exact interface without bogging you down in the Rust implementation details, here are the Python interface stub files that define the API.
> **`mbc_mom/geometry.pyi` (The Mesh Architecture):**
> ```python
> from typing import List
> 
> class Node:
>     x: float
>     y: float
>     z: float
>     def __init__(self, x: float, y: float, z: float) -> None: ...
> 
> class Segment:
>     start_idx: int
>     end_idx: int
>     radius: float
>     def __init__(self, start_node: int, end_node: int, radius: float) -> None: ...
> 
> class Dipole:
>     seg1_idx: int
>     seg2_idx: int
>     junction_idx: int
>     is_monopole: bool
> 
> class Mesh:
>     nodes: List[Node]
>     segments: List[Segment]
>     dipoles: List[Dipole]
> 
>     def __init__(self) -> None: ...
>     def add_node(self, node: Node) -> int: ...
>     def add_segment(self, segment: Segment) -> int: ...
>     def set_pec_ground(self) -> None: ...
>     def set_real_ground(self, sigma: float, eps_r: float, use_sommerfeld: bool = False) -> None: ...
>     def build_dipoles(self) -> None: ...
>     def validate(self, freq_hz: float) -> List[str]: ...
> 
> ```
> 
> 
> **`mbc_mom/__init__.pyi` (The Physics Solvers):**
> ```python
> from typing import List, Tuple
> from . import geometry
> 
> def compute_impedance_matrix(mesh: geometry.Mesh, freq_hz: float) -> List[complex]: ...
> 
> def compute_far_field(mesh: geometry.Mesh, currents: List[complex], freq_hz: float, thetas: List[float], phis: List[float]) -> List[float]: ...
> 
> def compute_near_field(
>     mesh: geometry.Mesh, 
>     currents: List[complex], 
>     freq_hz: float, 
>     xs: List[float], 
>     ys: List[float], 
>     zs: List[float]
> ) -> Tuple[List[complex], List[complex], List[complex], List[complex], List[complex], List[complex]]: ...
> 
> def get_incident_eval_points(
>     mesh: geometry.Mesh, 
>     points_per_seg: int = 7
> ) -> Tuple[List[float], List[float], List[float]]: ...
> 
> def compute_incident_v_matrix(
>     mesh: geometry.Mesh,
>     freq_hz: float,
>     ex: List[complex],
>     ey: List[complex],
>     ez: List[complex],
>     points_per_seg: int = 7
> ) -> List[complex]: ...
> 
> ```
> 
> 
> The engine features Multiradius Bridge-Currents (MBC), analytical near-field integration, Reflection Coefficient Approximation (RCA), exact Sommerfeld-Norton surface wave look-up tables, and incident field scattering capabilities. Please act as a senior technical writer and Python/Rust engineer to help me document this.
