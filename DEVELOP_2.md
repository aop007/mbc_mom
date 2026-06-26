Adding a ground plane is the natural next evolution of a professional Method of Moments solver. Transitioning from free-space to a half-space environment fundamentally changes the physics, transforming isolated wires into interacting arrays via image theory.

Integrating this into your hybrid Rust/Python architecture requires careful handling of both the topology (allowing current to flow *into* the ground) and the physics kernel (evaluating reflection coefficients and highly oscillatory integrals).

Here is the comprehensive architectural plan to implement Ground Planes in the MBC solver.

---

### 🗺️ Phase 1: Geometry and Data Structures (`src/geometry.rs`)

We must first define the physical properties of the ground and attach it to the `Mesh`.

**1. Define the Ground Enums and Structs:**

* Create a `GroundMaterial` enum: `PEC` or `Dielectric { sigma: f64, eps_r: f64 }`.
* Create a `GroundMethod` enum: `RCA` (Reflection Coefficient Approximation) or `Sommerfeld` (Exact).
* Create the `GroundPlane` struct holding these properties.

**2. Update the Mesh:**

* Add an `Option<GroundPlane>` field to the `Mesh` struct.
* Update the Python `#[pymethods]` to allow users to set the ground plane from the Python API (e.g., `mesh.set_pec_ground()` or `mesh.set_real_ground(sigma, eps_r, method)`).

---

### 🔗 Phase 2: Topology and the "Half-Dipole" Problem

In free space, current goes to zero at the physical ends of a wire. When a wire touches a ground plane (nodes at z=0), the current does *not* go to zero; it flows into the ground, mirrored by its image.

**1. Modify `build_dipoles` for Grounded Nodes:**

* Currently, a `Dipole` requires two physical segments.
* We must introduce a "Monopole Dipole" topology: If a segment ends at a node where z=0, the solver must treat the ground plane as the "second segment".
* The current basis function will peak at z=0 and smoothly transition into its virtual image below the ground.

---

### 🧮 Phase 3: The Impedance Matrix Kernel (`src/impedance.rs`)

This is where the physics gets heavy. Every mutual impedance calculation $Z_{ij}$ now consists of the direct interaction plus the interaction from the image of segment $j$.

**1. The Image Segment Generator:**

* Create a fast helper function that takes any `Segment` and generates its virtual image by flipping the Z-coordinates of its nodes and reversing the vertical vector components.

**2. Modify the Impedance Summation:**

* Update the nested parallel loops. Instead of just $Z_{ij} = Z_{direct}$, it becomes:

$$Z_{ij} = Z_{direct} + \Gamma \cdot Z_{image}$$


* **For PEC:** $\Gamma$ is strictly $+1$ for vertical currents and $-1$ for horizontal currents. The existing MBC integrals can be reused directly for the image segments.

**3. Implement Homogeneous Dielectric Approximations:**

* **Reflection Coefficient Approximation (RCA):** Implement the Fresnel reflection coefficients ($\Gamma_v$ and $\Gamma_h$) calculated from the wave polarization, $\epsilon_r$, and $\sigma$ evaluated at the target frequency. This is computationally cheap and works well for antennas slightly elevated above the ground.
* **Rigorous Sommerfeld Integrals:** For wires close to or touching lossy ground, RCA breaks down due to surface waves. Implement the exact Sommerfeld integrals. Because these are highly oscillatory and notoriously slow to compute, you will likely need to implement an interpolation lookup table (LUT) generated in Rust before the main matrix assembly begins.

---

### ⚡ Phase 4: Excitation and the Python API (`src/lib.rs` & Python)

The way we calculate input impedance changes when the feed is at the ground plane.

**1. Grounded Feeds (Monopoles):**

* When the Python script applies a $1.0\text{ V}$ Delta-Gap source to a junction at z=0, the engine is actually driving the gap between the physical segment and its image.
* The resulting input impedance $Z_{in}$ calculated by the Python script will be exactly half of the equivalent free-space dipole (e.g., a resonant quarter-wave monopole will correctly yield $\approx 36.5\ \Omega$ instead of $73\ \Omega$).

**2. Python Helper Functions:**

* Add specific feed-finding logic in Python to easily identify nodes located at $z \approx 0.0$ so users can easily drive grounded monopole antennas.

---

### 📡 Phase 5: Far-Field Adjustments (`src/farfield.rs`)

The radiation pattern of an antenna over a ground plane only exists in the upper hemisphere ($z > 0$).

**1. Superposition of Radiators:**

* During the precomputation of the radiating vector moments, automatically generate the image moments for every segment.
* Multiply the image moments by the appropriate $\Gamma$ (Fresnel coefficients for RCA, or $\pm 1$ for PEC) based on the observation angle $\theta$.

**2. Hemisphere Masking:**

* If the observation angle crosses below the horizon ($\theta > \pi/2$), force the computed radiation intensity $U(\theta, \phi)$ to exactly $0.0$.