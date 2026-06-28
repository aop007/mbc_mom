#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from mbc_mom import compute_impedance_matrix, compute_far_field, C
from mbc_mom.geometry import Mesh, Node, Segment


def build_dipole(mesh: Mesh, num_segments: int, length_m: float = 0.5, roll_deg: float = 0.0):
    roll = np.radians(roll_deg)
    
    # 0.5m Dipole oriented along the Z-axis
    z_coords = np.linspace(-length_m / 2, length_m / 2, num_segments + 1)
    nodes = [mesh.add_node(Node(float(z) * np.sin(roll), 0.0, float(z) * np.cos(roll))) for z in z_coords]
        
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.001))
    # end for
        
    return nodes[num_segments // 2]


def main():
    print("--- Computing Far-Field Radiation Pattern ---")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax = plt.subplot(111, polar=True)
    
    # Orient the plot so 0 degrees (Z-axis) is at the top
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    max_gain_list_dBi = []
    
    for dipole_length_m in [0.1, 0.25, 0.5, 0.75, 1.0]:
        mesh = Mesh()
        driven_node_idx = build_dipole(mesh, 30, length_m=dipole_length_m)
        mesh.build_dipoles()
        
        N = len(mesh.dipoles)
        freq_hz = C  # wavelength = 1 m
        
        # 1. Solve the Method of Moments
        Z = np.array(compute_impedance_matrix(mesh, freq_hz)).reshape((N, N))
        V = np.zeros(N, dtype=np.complex128)
        driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
        V[driven_dipole_idx] = 1.0 + 0.0j

        I = la.solve(Z, V)
        
        # 2. Setup Far-Field Grid
        # We want a full 360-degree Elevation cut (Theta from 0 to 2*pi)
        thetas = np.linspace(0, 2 * np.pi, 360).tolist()
        phis = [0.0] # Single azimuth slice for 2D plot
        
        # 3. Compute Radiation Intensity U(theta, phi) via Rust
        U_flat = compute_far_field(mesh, I.tolist(), freq_hz, thetas, phis)
        U = np.array(U_flat)
        
        # 4. Convert U to Directivity (dBi)
        # Total Radiated Power via input port
        Z_in = 1.0 / I[driven_dipole_idx]
        P_rad = 0.5 * np.real(V[driven_dipole_idx] * np.conj(I[driven_dipole_idx]))
        
        # D = 4 * pi * U / P_rad
        Directivity_linear = (4.0 * np.pi * U) / P_rad
        
        # Avoid log(0)
        Directivity_linear = np.maximum(Directivity_linear, 1e-10)
        Directivity_dBi = 10.0 * np.log10(Directivity_linear)
        
        max_dBi = np.max(Directivity_dBi)
        max_gain_list_dBi.append(max_dBi)
        
        print(f"Input Impedance: {Z_in.real:.2f} + j{Z_in.imag:.2f} Ohms")
        print(f"Max Directivity: {max_dBi:.2f} dBi @ {dipole_length_m:0.3f} m")
        ax.plot(thetas, Directivity_dBi, linewidth=2, label=f"{dipole_length_m:g} m")
    # end for

    # 5. Plotting

    abs_max_dBi = max(max_gain_list_dBi)

    # Clip the bottom of the plot at -30 dBi for clean visuals
    ax.set_ylim(-10, abs_max_dBi + 2)
    ax.grid()
    ax.legend()
    
    plt.title(f"Elevation Radiation Pattern (Gain: {abs_max_dBi:.2f} dBi)\nDipoles @ 300 MHz", va='bottom')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()