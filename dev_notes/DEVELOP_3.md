This is the Mount Everest of Computational Electromagnetics. Implementing rigorous Sommerfeld integrals will take your engine from a very good educational code to a true, industry-equivalent solver on par with the legendary NEC (Numerical Electromagnetics Code).

Because computing these integrals directly inside the $N \times N$ matrix loop would grind your CPU to a halt, we have to build a sophisticated "pre-processor" that maps out the physics before the matrix even begins assembling.

Here is our master plan to conquer the Sommerfeld Integrals.

---

### 🗺️ Phase 1: The Complex Integration Engine

We first need to build a standalone numerical integrator in Rust capable of handling the mathematical landmines inherent to Sommerfeld integrals.

* **The Integrand:** We will define the exact Sommerfeld correction term, which features a highly oscillatory Bessel function $J_0(\lambda \rho)$ and an exponential decay component.
* **Contour Deformation:** The denominator contains the "Sommerfeld Pole" (representing the surface wave). If we integrate straight along the real axis, the code will crash or diverge when it hits the pole. We must write a routine that mathematically deforms the integration path into the complex plane, loops *under* the pole, and comes back to the real axis.
* **Adaptive Quadrature:** We will implement a smart integrator that breaks the infinite upper bound ($0$ to $\infty$) into finite chunks between the zeros of the Bessel function, summing them until they converge.

### 📐 Phase 2: Asymptotic Bounds (The Safety Nets)

Numerical integration fails at extreme limits. We need fallback equations for the edges of our physics grid.

* **Near-Field Limit ($R \to 0$):** When the wire is microscopically close to the ground, the integrals blow up. We will implement Maclaurin series expansions to evaluate the fields here.
* **Far-Field Limit ($R \to \infty$):** At large distances, the integration takes too long to converge. We will implement **Norton's Asymptotic Expansions**, which perfectly approximate the ground wave over long distances algebraically.

### 🗄️ Phase 3: The 2D Lookup Table (LUT)

Once our integrator is bulletproof, we will build the pre-computation engine.

* **Grid Design:** We will define a 2D logarithmic grid based on radial distance ($\rho$) and the combined heights of the source and observation points ($z + z'$).
* **Parallel Generation:** Right before the MoM matrix starts, Rust will fire up all CPU cores to evaluate the Sommerfeld integrals at a few hundred key points on this grid, generating a map of the soil's exact electromagnetic response.
* **Bicubic Interpolation:** We will implement a fast 2D spline interpolator. When the matrix loop asks for a value between grid points, the LUT will instantly calculate a perfectly smooth estimation.

### ⚡ Phase 4: Splicing into the Matrix Kernel

We return to `src/impedance.rs` to replace the RCA approximation.

* Instead of computing a simple $\Gamma_v$ angle, the matrix will calculate $\rho$ and $z+z'$ for the interacting sub-segments.
* It will query the Sommerfeld LUT, instantly grabbing the exact surface wave coupling value.
* We add this value to the ideal image term, completing the exact Green's function for the lossy half-space.

---

**The Path Forward**
This is a complex operation, but we will take it one piece at a time. The most logical starting point is **Phase 1**: writing the raw mathematical function for the Sommerfeld integrand and plotting its poles in the complex plane to see what we are up against.

Shall we spin up a new Rust file specifically for the Sommerfeld math and start defining those complex contours?