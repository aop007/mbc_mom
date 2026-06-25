#!/usr/bin/env python


import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment


def build_half_wave_dipole(mesh: Mesh, num_segments: int, radius: float):
    """
    Builds a 0.5m length dipole for 300 MHz (lambda = 1m).
    """
    z_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = []
    
    for z in z_coords:
        nodes.append(mesh.add_node(Node(0.0, 0.0, float(z))))
        
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], radius))

    return nodes[num_segments // 2]

def generate_excitation_vector(mesh: Mesh, N: int, driven_node_idx: int, freq_hz: float, feed_type: str, a: float):
    """
    Generates the [V] excitation vector using either Delta-Gap or Magnetic Frill.
    """
    V = np.zeros(N, dtype=np.complex128)
    
    # 1. Delta-Gap Model
    if feed_type == "delta":
        driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
        V[driven_dipole_idx] = 1.0 + 0.0j
        return V, driven_dipole_idx

    # 2. Magnetic Frill Model
    # Assume a standard 50 Ohm coaxial feed geometry where b/a ~ 2.3
    b = a * 2.3 
    k = (2.0 * np.pi * freq_hz) / 299792458.0
    V_in = 1.0
    frill_factor = V_in / (2.0 * np.log(b / a))
    
    driven_node = mesh.nodes[driven_node_idx]
    z_feed = driven_node.z
    
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)

    # Numerically integrate the Frill E-field over every dipole in the mesh
    for i, dipole in enumerate(mesh.dipoles):
        seg1 = mesh.segments[dipole.seg1_idx]
        seg2 = mesh.segments[dipole.seg2_idx]
        
        junc_node = mesh.nodes[dipole.junction_idx]
        
        v_dipole = 0.0j
        
        # Integrate over both segments of the dipole (20-point midpoint rule)
        for seg, is_seg1 in [(seg1, True), (seg2, False)]:
            n_start = mesh.nodes[seg.start_idx]
            n_end = mesh.nodes[seg.end_idx]
            
            # Segment length and bounds
            L = np.sqrt((n_end.x - n_start.x)**2 + (n_end.y - n_start.y)**2 + (n_end.z - n_start.z)**2)
            z_min = min(n_start.z, n_end.z)
            z_max = max(n_start.z, n_end.z)
            
            steps = 20
            dz = L / steps
            
            for step in range(steps):
                z_eval = z_min + (step + 0.5) * dz
                
                # Distance from the feed point
                z_rel = z_eval - z_feed
                
                # Frill E-field at this point
                R_a = np.sqrt(z_rel**2 + a**2)
                R_b = np.sqrt(z_rel**2 + b**2)
                E_z = frill_factor * (np.exp(-1j * k * R_a) / R_a - np.exp(-1j * k * R_b) / R_b)
                
                # PWS Test Function weighting (Distance from the dipole's own junction)
                dist_to_junc = np.abs(z_eval - junc_node.z)
                T_z = np.sin(k * (L - dist_to_junc)) / np.sin(k * L)
                
                # Directional sign based on KCL continuous flow
                junc_is_end = (seg.end_idx == dipole.junction_idx)
                sign = 1.0 if (is_seg1 == junc_is_end) else -1.0
                
                v_dipole += sign * E_z * T_z * dz
                
        V[i] = v_dipole
        
    return V, driven_dipole_idx

def run_analysis(num_segments: int, feed_type: str):
    mesh = Mesh()
    radius = 0.001
    freq_hz = 300e6 
    
    driven_node_idx = build_half_wave_dipole(mesh, num_segments, radius)
    mesh.build_dipoles()
    N = len(mesh.dipoles)
    
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    V, driven_dipole_idx = generate_excitation_vector(mesh, N, driven_node_idx, freq_hz, feed_type, radius)

    I = la.solve(Z, V)
    
    I_in = I[driven_dipole_idx]
    
    # For the frill, V[driven_dipole] is the projected voltage, but the applied coax voltage is V_in = 1.0
    V_applied = 1.0 if feed_type == "frill" else V[driven_dipole_idx]
    Z_in = V_applied / I_in

    print(f"{num_segments:8d} | {Z_in.real:10.4f} | {Z_in.imag:10.4f}")

def main():
    print("--- Convergence Analysis: Delta vs Frill (R=1mm) ---")
    print(f"{'Segments':>8} | {'Re(Z_in)':>10} | {'Im(Z_in)':>10}")
    
    print("\n[ Delta-Gap Feed ]")
    for segs in [2, 10, 30, 60, 98]:
        run_analysis(segs, "delta")
        
    print("\n[ Magnetic Frill Feed ]")
    for segs in [2, 10, 30, 60, 98]:
        run_analysis(segs, "frill")

if __name__ == "__main__":
    main()