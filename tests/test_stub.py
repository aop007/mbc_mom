#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment

def build_stepped_dipole(mesh: Mesh, forward: bool):
    """
    Builds a 0.75m dipole centered at Z=0.
    Radius A = 1 mm
    Radius B = 50 mm (50:1 step ratio)
    """
    r_thin = 0.001
    r_thick = 0.050

    # Assign radii based on drawing direction to test reciprocity
    r1 = r_thin if forward else r_thick
    r2 = r_thick if forward else r_thin

    n0 = mesh.add_node(Node(0.0, 0.0, -0.375))
    n1 = mesh.add_node(Node(0.0, 0.0,  0.000)) # The Junction
    n2 = mesh.add_node(Node(0.0, 0.0,  0.375))

    mesh.add_segment(Segment(n0, n1, r1))
    mesh.add_segment(Segment(n1, n2, r2))

    return n1

def test_multiradius_reciprocity():
    freq_hz = 100e6  # 100 MHz
    
    print("--- MBC Multiradius Junction Test (50:1 Ratio) ---")
    
    for direction in [True, False]:
        mesh = Mesh()
        driven_node_idx = build_stepped_dipole(mesh, forward=direction)
        mesh.build_dipoles()
        
        N = len(mesh.dipoles)
        
        # Compute the Z matrix in Rust
        flat_z = compute_impedance_matrix(mesh, freq_hz)
        Z = np.array(flat_z).reshape((N, N))
        
        # Delta-Gap Excitation
        V = np.zeros(N, dtype=np.complex128)
        driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
        V[driven_dipole_idx] = 1.0 + 0.0j
        
        # Solve
        I = la.solve(Z, V)
        Z_in = 1.0 / I[driven_dipole_idx]
        
        # Display the geometry config and result
        config = "Thin -> Thick" if direction else "Thick -> Thin"
        print(f"\nGeometry: {config}")
        print(f"Z_11 Matrix Term:  {Z[0, 0].real:.6f} + {Z[0, 0].imag:.6f}j")
        print(f"Input Impedance:   {Z_in.real:.6f} + {Z_in.imag:.6f}j Ohms")

if __name__ == "__main__":
    test_multiradius_reciprocity()