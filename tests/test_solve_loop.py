#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment


def build_rectangular_loop(mesh: Mesh):
    """
    Builds a 30 mm x 7.5 mm loop with wire radius 1.25 mm.
    We split each side into 2 segments so we have center nodes to drive.
    """
    r = 0.00125  # 1.25 mm radius
    
    # 1. Define Nodes (converted to meters)
    nodes = [
        mesh.add_node(Node(0.0,    0.0,     0.0)), # 0: Bottom-left
        mesh.add_node(Node(0.015,  0.0,     0.0)), # 1: Bottom-center (Driven)
        mesh.add_node(Node(0.030,  0.0,     0.0)), # 2: Bottom-right
        mesh.add_node(Node(0.030,  0.00375, 0.0)), # 3: Right-center
        mesh.add_node(Node(0.030,  0.0075,  0.0)), # 4: Top-right
        mesh.add_node(Node(0.015,  0.0075,  0.0)), # 5: Top-center
        mesh.add_node(Node(0.00,   0.0075,  0.0)), # 6: Top-left
        mesh.add_node(Node(0.00,   0.00375, 0.0)), # 7: Left-center
    ]

    # 2. Define Segments forming the closed loop
    mesh.add_segment(Segment(nodes[0], nodes[1], r))
    mesh.add_segment(Segment(nodes[1], nodes[2], r))
    mesh.add_segment(Segment(nodes[2], nodes[3], r))
    mesh.add_segment(Segment(nodes[3], nodes[4], r))
    mesh.add_segment(Segment(nodes[4], nodes[5], r))
    mesh.add_segment(Segment(nodes[5], nodes[6], r))
    mesh.add_segment(Segment(nodes[6], nodes[7], r))
    mesh.add_segment(Segment(nodes[7], nodes[0], r))

    return nodes[1]  # Return the index of our driven junction

def main():
    print("--- MBC Method of Moments Solver ---")
    mesh = Mesh()
    
    # 1. Build Geometry
    driven_node_idx = build_rectangular_loop(mesh)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    print(f"Generated Mesh: {len(mesh.segments)} segments, {N} dipoles.")

    # 2. Generate Impedance Matrix [Z] (100 MHz)
    freq_hz = 100e6
    print(f"Computing Z matrix at {freq_hz / 1e6} MHz in Rust (Parallel)...")
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    # 3. Generate Excitation Vector [V]
    # We apply a 1.0V Delta-Gap source at the dipole located at 'driven_node_idx'
    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = None
    
    for i, dipole in enumerate(mesh.dipoles):
        if dipole.junction_idx == driven_node_idx:
            driven_dipole_idx = i
            break
            
    if driven_dipole_idx is None:
        raise ValueError("Could not find a dipole at the driven junction!")
        
    V[driven_dipole_idx] = 1.0 + 0.0j

    # 4. Solve the Linear System: [Z][I] = [V]
    print("Solving linear system [Z][I] = [V]...")
    I = la.solve(Z, V)

    # 5. Extract Input Impedance
    I_in = I[driven_dipole_idx]
    Z_in = 1.0 / I_in

    print("\n--- Results ---")
    print(f"Driven Dipole Index: {driven_dipole_idx}")
    print(f"Input Current (I):   {I_in:.4e} A")
    print(f"Input Impedance (Z): {Z_in.real:.6e} + {Z_in.imag:.4f}j Ohms")

if __name__ == "__main__":
    main()