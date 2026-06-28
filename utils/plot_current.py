#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix, C
from mbc_mom.geometry import Mesh, Node, Segment

def build_dipole(mesh: Mesh, num_segments: int):
    """Builds a 0.5m half-wave dipole centered at Z=0."""
    z_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = []
    
    for z in z_coords:
        nodes.append(mesh.add_node(Node(0.0, 0.0, float(z))))
        
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.001))
        
    return nodes[num_segments // 2]

def main():
    print("--- Extracting Standing Wave Current Distribution ---")
    mesh = Mesh()
    num_segments = 30
    
    driven_node_idx = build_dipole(mesh, num_segments)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = C  # wavelength = 1 m
    
    # 1. Compute and Solve
    flat_z = compute_impedance_matrix(mesh, freq_hz)
    Z = np.array(flat_z).reshape((N, N))

    V = np.zeros(N, dtype=np.complex128)
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    I = la.solve(Z, V)

    # 2. Extract Data for Plotting
    i_mags_z_m = []
    
    # The solved I vector holds the current exactly at the junction of each dipole
    for idx, dipole in enumerate(mesh.dipoles):
        junc_z = mesh.nodes[dipole.junction_idx].z
        
        # Convert to milliamps for better readability
        i_mags_z_m.append((np.abs(I[idx]) * 1000.0, junc_z)) 
    # end for
    
    i_mags_z_m = list(sorted(i_mags_z_m, key=lambda x: x[1]))

    z_vals = [z for _, z in i_mags_z_m]
    i_mags = [i_mA for i_mA, _ in i_mags_z_m]

    # 3. Enforce Boundary Conditions at the wire tips
    # Current must go to zero at the physical ends of the dipole
    z_vals = [-0.25] + z_vals + [0.25]
    i_mags = [0.0] + i_mags + [0.0]

    # 4. Plot
    plt.figure(figsize=(10, 5))
    plt.plot(z_vals, i_mags, marker='o', linestyle='-', color='#1f77b4', linewidth=2, markersize=5)
    plt.fill_between(z_vals, i_mags, alpha=0.2, color='#1f77b4')
    
    plt.title(f"Standing Wave on Half-Wave Dipole (300 MHz, {num_segments} segments)", fontsize=14)
    plt.xlabel("Position along Wire (meters)", fontsize=12)
    plt.ylabel("Current Magnitude |I| (mA)", fontsize=12)
    
    plt.axvline(0, color='red', linestyle='--', label='Driven Feed Point')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.xlim(-0.26, 0.26)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()