#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment

def build_monopole(mesh: Mesh, num_segments: int):
    """
    Builds a 0.25m monopole (quarter-wave at 300 MHz) starting at Z=0.
    """
    radius = 0.0001  # Thin wire to match analytical King-Middleton approximations
    
    # Create nodes along the Z axis from 0.0 up to 0.25
    z_coords = np.linspace(0.0, 0.25, num_segments + 1)
    nodes = []
    
    for z in z_coords:
        nodes.append(mesh.add_node(Node(0.0, 0.0, float(z))))
        
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], radius))

    # Return the index of the base node (touching the ground plane)
    return nodes[0]

def main():
    print("--- MBC Solver: Quarter-Wave Monopole over PEC Ground ---")
    mesh = Mesh()
    
    # 1. Activate the PEC Ground Plane BEFORE building dipoles
    mesh.set_pec_ground()
    
    # 2. Build Geometry
    num_segments = 15
    driven_node_idx = build_monopole(mesh, num_segments)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    print(f"Generated Mesh: {len(mesh.segments)} physical segments, {N} dipoles.")
    
    # Check if the solver correctly identified the monopole
    driven_dipole = next(d for d in mesh.dipoles if d.junction_idx == driven_node_idx)
    print(f"Base Dipole identified as grounded monopole: {driven_dipole.is_monopole}")

    # 3. Compute and Solve
    freq_hz = 300e6  # 300 MHz
    print(f"Computing Z matrix at {freq_hz / 1e6} MHz (including virtual images)...")
    
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    # Drive the gap between the physical wire and its virtual image
    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    print("Solving linear system [Z][I] = [V]...")
    I = la.solve(Z, V)

    I_in = I[driven_dipole_idx]
    Z_in = 1.0 / I_in

    print("\n--- Results ---")
    print(f"Input Impedance (Z): {Z_in.real:.2f} + {Z_in.imag:.2f}j Ohms")
    print(f"Expected Target:     ~36.55 + 21.25j Ohms (Half of a free-space dipole)")

if __name__ == "__main__":
    main()