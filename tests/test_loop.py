#!/usr/bin/env python

import matplotlib.pyplot as pyplot
import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix, C
from mbc_mom.geometry import Mesh, Node, Segment

def build_rectangular_loop(mesh: Mesh):
    """
    Builds a 30 mm x 7.5 mm loop with a wire radius of 1.25 mm.
    """
    r = 0.00125  
    
    nodes = [
        mesh.add_node(Node(0.0,    0.0,     0.0)),     # 0: Bottom-left
        mesh.add_node(Node(0.015,  0.0,     0.0)),     # 1: Bottom-center (Driven)
        mesh.add_node(Node(0.030,  0.0,     0.0)),     # 2: Bottom-right
        mesh.add_node(Node(0.030,  0.00375, 0.0)),     # 3: Right-center
        mesh.add_node(Node(0.030,  0.0075,  0.0)),     # 4: Top-right
        mesh.add_node(Node(0.015,  0.0075,  0.0)),     # 5: Top-center
        mesh.add_node(Node(0.00,   0.0075,  0.0)),     # 6: Top-left
        mesh.add_node(Node(0.00,   0.00375, 0.0)),     # 7: Left-center
    ]

    # Close the loop
    mesh.add_segment(Segment(nodes[0], nodes[1], r))
    mesh.add_segment(Segment(nodes[1], nodes[2], r))
    mesh.add_segment(Segment(nodes[2], nodes[3], r))
    mesh.add_segment(Segment(nodes[3], nodes[4], r))
    mesh.add_segment(Segment(nodes[4], nodes[5], r))
    mesh.add_segment(Segment(nodes[5], nodes[6], r))
    mesh.add_segment(Segment(nodes[6], nodes[7], r))
    mesh.add_segment(Segment(nodes[7], nodes[0], r))

    return nodes[1] 

def main():
    print("--- MBC Solver: Small Rectangular Loop Benchmark ---")
    mesh = Mesh()
    driven_node_idx = build_rectangular_loop(mesh)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = 100e6  # 100 MHz
    
    print(f"Mesh: {len(mesh.segments)} segments, {N} dipoles.")
    print(f"Computing Z matrix at {freq_hz / 1e6} MHz...")
    
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    # Delta-gap feed
    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    I = la.solve(Z, V)
    Z_in = 1.0 / I[driven_dipole_idx]

    # Extract Physics Metrics
    r_rad_micro = Z_in.real * 1e6
    omega = 2.0 * np.pi * freq_hz
    inductance_nh = (Z_in.imag / omega) * 1e9

    print("\n--- Results ---")
    print(f"Input Impedance:      {Z_in.real:.8e} + {Z_in.imag:.4f}j Ohms")
    print(f"Radiation Resistance: {r_rad_micro:.2f} micro-Ohms")
    print(f"Equivalent Inductance:{inductance_nh:.2f} nH")
    
    # Analytical verification
    area = 0.030 * 0.0075
    wavelength = C / freq_hz
    r_analytical = 31200.0 * (area / (wavelength**2))**2 * 1e6
    print(f"\nAnalytical Target:    ~{r_analytical:.2f} micro-Ohms")
    
    fig, ax = pyplot.subplots()
    
    for segment in mesh.segments:
        start_node = mesh.nodes[segment.start_idx]
        stop_node = mesh.nodes[segment.end_idx]
        
        ax.plot(
            [start_node.x, stop_node.x],
            [start_node.y, stop_node.y],
            marker='o',
        )
    # end for
    
    ax.grid()
    ax.legend()
    
    pyplot.show()

if __name__ == "__main__":
    main()