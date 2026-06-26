#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from mbc_mom import compute_impedance_matrix, compute_far_field
from mbc_mom.geometry import Mesh, Node, Segment

def build_monopole(mesh: Mesh, num_segments: int, height_m: float = 0.25):
    z_coords = np.linspace(0.0, height_m, num_segments + 1)
    nodes = [mesh.add_node(Node(0.0, 0.0, float(z))) for z in z_coords]
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.0001))
    return nodes[0]

def main(use_pec: bool = False, use_dry_gnd: bool = False, use_saltwater: bool = False):
    mesh = Mesh()

    if use_pec:
        print("--- Far-Field Pattern: Monopole over PEC Ground ---")
        mesh.set_pec_ground() # Turn on the ground!
    elif use_dry_gnd:
        print("--- Far-Field Pattern: Monopole over Dry Ground Water ---")
        mesh.set_real_ground(sigma=0.001, eps_r=3.0, use_sommerfeld=False)
    elif use_saltwater:
        print("--- Far-Field Pattern: Monopole over Salt Ground ---")
        mesh.set_real_ground(sigma=5.0, eps_r=81.0, use_sommerfeld=False)
    else:
        print("--- Far-Field Pattern: Monopole over Void ---")
    # end if
    
    driven_node_idx = build_monopole(mesh, 15, height_m=0.25)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = 300e6
    
    # 1. Solve MoM
    Z = np.array(compute_impedance_matrix(mesh, freq_hz)).reshape((N, N))
    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    I = la.solve(Z, V)
    
    # 2. Far-Field Grid (Elevation cut from 0 to 180 degrees)
    # Even though we calculate to 180, the Rust kernel will mask > 90 to zero
    thetas = np.linspace(0, np.pi, 180).tolist()
    phis = [0.0]
    
    U_flat = compute_far_field(mesh, I.tolist(), freq_hz, thetas, phis)
    U = np.array(U_flat)
    
    # 3. Directivity Calculation
    # For a grounded monopole, input power is half that of a dipole
    Z_in = 1.0 / I[driven_dipole_idx]
    P_rad = 0.5 * np.real(V[driven_dipole_idx] * np.conj(I[driven_dipole_idx]))
    
    Directivity_linear = (4.0 * np.pi * U) / P_rad
    Directivity_linear = np.maximum(Directivity_linear, 1e-10)
    Directivity_dBi = 10.0 * np.log10(Directivity_linear)
    
    max_dBi = np.max(Directivity_dBi)
    print(f"Input Impedance: {Z_in.real:.2f} + {Z_in.imag:.2f}j Ohms")
    print(f"Max Directivity: {max_dBi:.2f} dBi")

    # 4. Plot
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    # We mirror the positive theta values to the left side of the plot for a full slice
    thetas_full = np.concatenate([-np.array(thetas)[::-1], thetas])
    D_full = np.concatenate([Directivity_dBi[::-1], Directivity_dBi])
    
    ax.plot(thetas_full, D_full, color='blue', linewidth=2)
    ax.set_ylim(-30, max_dBi + 2)
    
    # Draw a line representing the ground plane
    ax.axhline(0, color='black', linewidth=3)
    
    plt.title(f"Grounded Quarter-Wave Monopole\nGain: {max_dBi:.2f} dBi", va='bottom')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main(use_pec=True)
    main(use_dry_gnd=True)
    main(use_saltwater=True)