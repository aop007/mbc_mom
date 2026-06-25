# mbc-mom: Multiradius Bridge-Current Method of Moments

A high-performance computational electromagnetics (CEM) solver implementing the Multiradius Bridge-Current (MBC) Method of Moments. 

This solver builds upon the foundational piecewise-sinusoidal thin-wire codes developed by J.H. Richmond and extends them using the formulation by Mark A. Tilston and Keith G. Balmain. The MBC method features an exactly symmetric mutual impedance matrix, ensuring reciprocity between sources. It is completely unconstrained regarding both the length ratio and the radius ratio of adjoining segments, provided the wires remain electrically thin.

---

## 🏗️ Project Architecture

This is a **Mixed Python/Rust** project. It combines the high-level ease of Python for geometry definition and post-processing with the raw performance of Rust for dense complex matrix filling and mathematical integrations.

* **`src/lib.rs` (`mbc_lib`)**: The compiled Rust backend containing the memory-safe data structures (`Node`, `Segment`, `Mesh`) and the parallelized impedance matrix solvers. Compiled via PyO3.
* **`mbc_mom/`**: The frontend Python package. It acts as a wrapper around `mbc_lib`, exposing a clean API and providing static typing signatures (`.pyi` files) for language servers like Pylance/Jedi.
* **`pyproject.toml`**: Configured to use `maturin` as the build backend, injecting the compiled Rust `.so` directly into the `mbc_mom` Python package.

---

## 🚀 Environment & Setup

This project uses `pixi` for hermetic, cross-platform dependency management, completely isolating the Rust toolchain, Python, and `maturin` from the host OS (e.g., WSL Ubuntu).

### 1. Initializing the Environment
If cloning this repository fresh, install the dependencies lockfile:
```bash
pixi install
```

### 2. Building the Rust Extension

```bash
pixi run build-dev
```

### 3. Cleaning the environment

```bash
pixi run pip uninstall mbc-mom -y
pixi run cargo clean
pixi run build-dev
```

### 4. Running tests

```bash
pixi run python -m tests.test_rust_interface
```

