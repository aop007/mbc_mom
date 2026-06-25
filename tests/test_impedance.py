#!/usr/bin/env python

import tabulate
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as pyplot

from mbc_mom.geometry import Mesh, Node, Segment
from mbc_mom import compute_impedance_matrix


def test_parallel_dipoles():
    mesh = Mesh()

    # Step-radius junction
    n0 = mesh.add_node(Node(0.0, 0.0, -0.5))
    n1 = mesh.add_node(Node(0.0, 0.0, 0.0))  # Junction node
    n2 = mesh.add_node(Node(0.0, 0.0, 0.5))

    mesh.add_segment(Segment(n0, n1, radius=0.001))
    mesh.add_segment(Segment(n1, n2, radius=0.005))

    # Add parallel VDipole
    n3 = mesh.add_node(Node(1.0, 0.0, -0.5))
    n4 = mesh.add_node(Node(1.0, 0.0, 0.0))  # Junction node
    n5 = mesh.add_node(Node(1.0, 0.0, 0.5))

    mesh.add_segment(Segment(n3, n4, radius=0.001))
    mesh.add_segment(Segment(n4, n5, radius=0.005))

    mesh.build_dipoles()

    # Get the flat matrix from Rust
    flat_z = compute_impedance_matrix(mesh, 100e6) # 100 MHz

    # Reshape into a 2D NxN NumPy array
    N = len(mesh.dipoles)
    Z_matrix = np.array(flat_z).reshape((N, N))

    print(f"Matrix shape: {Z_matrix.shape}")
    print("Exact Reciprocity Check (Z_01 == Z_10):", Z_matrix[0,1] == Z_matrix[1,0])
# end test_parallel_dipoles()

def build_half_wave_dipole(mesh: Mesh, num_segments: int) -> int:
    """
    Builds a 0.5m length dipole for 300 MHz (lambda = 1m).
    """
    r = 0.001  # Very thin wire to approximate King-Middleton
    
    # Create nodes along the Z axis from -0.25 to 0.25
    z_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = []
    
    for z in z_coords:
        nodes.append(mesh.add_node(Node(0.0, 0.0, float(z))))
        
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], r))
    # end for

    # Return the index of the center node (the feed point)
    return nodes[num_segments // 2]
# end build_half_wave_dipole()

def test_half_wave_dipole(
    num_segments: int = 2,
) -> complex:
    print("--- MBC Solver: Half-Wave Dipole ---")
    
    expected_input_Z_Ohms = 73.1 + 42.5j
    
    mesh = Mesh()
    driven_node_idx = build_half_wave_dipole(mesh=mesh, num_segments=num_segments)
    
    driven_node = mesh.nodes[driven_node_idx]
    assert np.abs(driven_node.z) < 0.01, f"driven_node.z: {driven_node.z}"
    
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    print(f"Generated Mesh: {len(mesh.segments)} segments, {N} dipoles.")
    
    assert N % 2 == 1, f"Odd number of dipoles is expected!"

    # 300 MHz
    freq_hz = 300e6 
    print(f"Computing Z matrix at {freq_hz / 1e6} MHz...")
    
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    excitation_Vp = 1.0 + 0.0j 
    V[driven_dipole_idx] = excitation_Vp

    print("Solving linear system [Z][I] = [V]...")
    
    try:
        I = la.solve(Z, V)
    except Exception:
        print(f"Z: {Z}")
        print(f"V: {V}")
        raise
    # end try

    I_in = I[driven_dipole_idx]
    Z_in = excitation_Vp / I_in
    
    error_Z_Ohms = Z_in - expected_input_Z_Ohms

    print(f"\n--- Results {num_segments} segments ---")
    print(f"Input Impedance (Z): {Z_in.real:.2f} + {Z_in.imag:.2f}j Ohms")
    print(f"Expected (King-Middleton): ~73.10 + 42.50j Ohms Error: {error_Z_Ohms.real:.2f} + {error_Z_Ohms.imag:.2f}j Ohms")
    
    return Z_in
# end test_half_wave_dipole()

if __name__ == "__main__":
    num_segment_array = np.arange(2, 102, 4)
    
    z_in_list = np.zeros_like(num_segment_array, dtype=complex)
    
    for ix, num_segments in enumerate(num_segment_array):
        Z_in = test_half_wave_dipole(num_segments=int(num_segments))
        
        z_in_list[ix] = Z_in
    # end for
    
    print(
        tabulate.tabulate(
            [[segs, np.real(Z_in), np.imag(Z_in)] for segs, Z_in in zip(num_segment_array, z_in_list)],
            headers="Segments|Re(Z_in)|Im(Z_in)|".split('|')
        )
    )
    
    fig, ax = pyplot.subplots()
    
    ax.plot(
        num_segment_array,
        np.real(z_in_list),
        marker='o',
        label='Re(Z_in)',
    )
    
    ax.plot(
        num_segment_array,
        np.imag(z_in_list),
        marker='x',
        label='Ie(Z_in)',
    )
    
    ax.grid()
    ax.legend()
    ax.semilogx()
    
    pyplot.show()
# end if