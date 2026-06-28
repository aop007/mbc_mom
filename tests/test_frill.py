#!/usr/bin/env python


import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from mbc_mom import compute_impedance_matrix, C
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
    k = (2.0 * np.pi * freq_hz) / C
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

def run_analysis(num_segments: int, feed_type: str) -> complex:
    mesh = Mesh()
    radius = 0.001
    freq_hz = C  # wavelength = 1 m 
    
    driven_node_idx = build_half_wave_dipole(mesh, num_segments, radius)
    mesh.build_dipoles()
    N = len(mesh.dipoles)
    
    warnings = mesh.validate(freq_hz)
    if warnings:
        print("⚠️ Mesh Validation Warnings:")
        for w in warnings:
            print(f"  - {w}")
        # end for
        print("⚠️ *************************")
        
        raise RuntimeError("Found warnings!")
    # end if
    
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    V, driven_dipole_idx = generate_excitation_vector(mesh, N, driven_node_idx, freq_hz, feed_type, radius)

    I = la.solve(Z, V)
    
    I_in = I[driven_dipole_idx]
    
    # For the frill, V[driven_dipole] is the projected voltage, but the applied coax voltage is V_in = 1.0
    V_applied = 1.0 if feed_type == "frill" else V[driven_dipole_idx]
    Z_in = V_applied / I_in

    print(f"{num_segments:8d} | {Z_in.real:10.4f} | {Z_in.imag:10.4f}")
    
    return Z_in

def main(show=False):
    print("--- Convergence Analysis: Delta vs Frill (R=1mm) ---")
    print(f"{'Segments':>8} | {'Re(Z_in)':>10} | {'Im(Z_in)':>10}")
    
    segments = 2.0 ** np.arange(1, 8, 0.125)
    segments_int = np.asarray(segments, dtype=int)
    segments_int = segments_int[np.where(segments_int % 2 == 0)]  # Keep even nb of segments
    segments_int = np.array(list(sorted(set(segments_int.tolist()))))
    
    print("\n[ Delta-Gap Feed ]")
    
    delta_gap_z_in = [run_analysis(segs, "delta") for segs in segments_int]
        
    print("\n[ Magnetic Frill Feed ]")
    
    mag_frill_gap_z_in = [run_analysis(segs, "frill") for segs in segments_int]
        
    fig, ax = plt.subplots()
    
    for op, linestyle, prefix in [
        (np.real, '-', 'Re'),
        (np.imag, '--', 'Im'),
    ]:
        ax.plot(
            segments_int,
            op(delta_gap_z_in),
            marker='o',
            linestyle=linestyle,
            label=f'{prefix}(Delta Gap)',
        )
        
        ax.plot(
            segments_int,
            op(mag_frill_gap_z_in),
            marker='o',
            linestyle=linestyle,
            label=f'{prefix}(Magnetic Frill)',
        )
    # end for
    
    ax.grid()
    ax.legend()
    ax.semilogx()
    ax.set_title("Input Impedance vs Segments")
    ax.set_xlabel("Segments")
    ax.set_ylabel("$Z_{in}$")
    
    if show:
        plt.show()
    # end if
# end main()

if __name__ == "__main__":
    main(show=True)