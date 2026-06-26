#!/usr/bin/env python

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mbc_mom import compute_impedance_matrix, compute_near_field
from mbc_mom.geometry import Mesh, Node, Segment

def build_dipole(mesh: Mesh, num_segments: int):
    # 0.5m Half-Wave Dipole at 300 MHz, oriented along the Z-axis
    z_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = [mesh.add_node(Node(0.0, 0.0, float(z))) for z in z_coords]
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.001))
    return nodes[num_segments // 2]

def main():
    print("--- Computing Reactive Near-Field ---")
    mesh = Mesh()
    driven_node_idx = build_dipole(mesh, 30)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = 300e6
    
    # 1. Solve MoM for the exact current distribution
    print("Solving Impedance Matrix...")
    Z = np.array(compute_impedance_matrix(mesh, freq_hz)).reshape((N, N))
    V = np.zeros(N, dtype=np.complex128)
    
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    I = la.solve(Z, V)
    
    # 2. Setup a 2D Spatial Grid (X-Z plane, slicing right through the antenna)
    print("Generating 10,000-point Spatial Grid...")
    grid_res = 100
    x_vals = np.linspace(-0.5, 0.5, grid_res)
    z_vals = np.linspace(-0.5, 0.5, grid_res)
    
    X, Z_grid = np.meshgrid(x_vals, z_vals)
    Y = np.zeros_like(X) # Y = 0 slice
    
    # Flatten grids to 1D lists to eliminate Python loop overhead
    xs_flat = X.flatten().tolist()
    ys_flat = Y.flatten().tolist()
    zs_flat = Z_grid.flatten().tolist()
    
    # 3. Call the parallelized Rust kernel
    print("Calculating Exact E and H fields in Rust...")
    Ex_flat, Ey_flat, Ez_flat, Hx_flat, Hy_flat, Hz_flat = compute_near_field(
        mesh, I.tolist(), freq_hz, xs_flat, ys_flat, zs_flat
    )
    
    # 4. Reshape back to 2D NumPy arrays
    Ex = np.array(Ex_flat).reshape(grid_res, grid_res)
    Ey = np.array(Ey_flat).reshape(grid_res, grid_res)
    Ez = np.array(Ez_flat).reshape(grid_res, grid_res)
    
    # Calculate Total Electric Field Magnitude: |E| = sqrt(|Ex|^2 + |Ey|^2 + |Ez|^2)
    E_mag = np.sqrt(np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2)
    
    Hx = np.array(Hx_flat).reshape(grid_res, grid_res)
    Hy = np.array(Hy_flat).reshape(grid_res, grid_res)
    Hz = np.array(Hz_flat).reshape(grid_res, grid_res)
    
    # Calculate Total Electric Field Magnitude: |E| = sqrt(|Ex|^2 + |Ey|^2 + |Ez|^2)
    H_mag = np.sqrt(np.abs(Hx)**2 + np.abs(Hy)**2 + np.abs(Hz)**2)
    
    # 5. Plot the fields using Matplotlib
    print("Plotting results...")
    fig, [ax_E, ax_H] = plt.subplots(figsize=(8, 7), ncols=2)
    
    for ax, field, title in [(ax_E, E_mag, 'Electric Field Magnitude |E| (V/m)'), (ax_H, H_mag, 'Magnetic Field Magnitude |H| (A/m)')]:
        # LogNorm handles the massive dynamic range (from V/m to kV/m) cleanly
        vmax=np.max(field)
        
        pcm = ax.pcolormesh(
            X, 
            Z_grid, 
            field, 
            shading='auto', 
            cmap='inferno', 
            norm=LogNorm(vmin=vmax / 1000.0, vmax=vmax)
        )
        plt.colorbar(pcm, label=title)
        
        # Draw the physical antenna wire as a thick white line
        ax.plot([0, 0], [-0.25, 0.25], color='white', linewidth=3, solid_capstyle='round', label='Dipole Antenna')
    
        ax.set_title('Near-Field Reactive Energy (|E|/|H|)\n300 MHz Half-Wave Dipole')
        ax.set_xlabel('X (meters)')
        ax.set_ylabel('Z (meters)')
        ax.axis('equal') # Ensure true spatial proportions
        ax.legend(loc='upper right')
    # end for
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()