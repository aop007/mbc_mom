#!/usr/bin/env python

import matplotlib.pyplot as pyplot
import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment

def build_resonant_stub(mesh: Mesh, a1: float, a2: float):
    top_nodes = []
    bot_nodes = []
    
    # Length 750 mm, centered at x=0
    x_coords = np.linspace(-0.375, 0.375, 21)
    
    for x in x_coords:
        # Width 7.5 mm in the z-axis (top at +3.75mm, bottom at -3.75mm)
        top_nodes.append(mesh.add_node(Node(float(x), 0.0, 0.00375)))
        bot_nodes.append(mesh.add_node(Node(float(x), 0.0, -0.00375)))
        
    for i in range(20):
        mesh.add_segment(Segment(top_nodes[i], top_nodes[i+1], a1))
        mesh.add_segment(Segment(bot_nodes[i], bot_nodes[i+1], a2))
        
    # Shorting segment at left end (x < 0) using the nominal 1.25 mm radius
    mesh.add_segment(Segment(top_nodes[0], bot_nodes[0], 0.00125))
    
    # Return the indices for the junctions at x=0
    return top_nodes[10], bot_nodes[10]


def main():
    print("--- MBC Solver: Resonant Two-Wire Stub Benchmark ---")
    print(f"{'log10(A1/A2)':>12} | {'A1 (mm)':>10} | {'A2 (mm)':>10} | {'Re(Z_in)':>12} | {'Im(Z_in)':>12}")
    print("-" * 65)
    
    freq_hz = 99.93e6
    
    for log_r in range(-6, 7):
        # Apply the specific A1/A2 ratio rules
        if log_r < 0:
            a2 = 0.00125
            a1 = a2 * (10.0 ** log_r)
        else:
            a1 = 0.00125
            a2 = a1 / (10.0 ** log_r)
        # end if
        
        mesh = Mesh()
        driven_top_idx, driven_bot_idx = build_resonant_stub(mesh, a1, a2)
        mesh.build_dipoles()
        
        N = len(mesh.dipoles)
        flat_z = compute_impedance_matrix(mesh, freq_hz)
        Z = np.array(flat_z).reshape((N, N))
        
        V = np.zeros(N, dtype=np.complex128)
        
        # Locate the dipoles spanning the feed nodes
        top_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_top_idx)
        bot_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_bot_idx)
        
        # KCL Feed Setup:
        # Top segment positive end towards x > 0 (Applied Voltage = +1.0V)
        # Bottom segment positive end towards x < 0 (Applied Voltage = -1.0V)
        V[top_dipole_idx] = 1.0 + 0.0j
        V[bot_dipole_idx] = -1.0 + 0.0j
        
        I = la.solve(Z, V)
        
        # Extract currents
        I_top = I[top_dipole_idx]
        I_bot = I[bot_dipole_idx]
        
        # Differential Input Impedance
        # V_diff = V_top - V_bot = 2.0V
        # I_diff is the average loop current
        I_diff = (I_top - I_bot) / 2.0
        Z_in = 2.0 / I_diff
        
        print(f"{log_r:12d} | {a1*1000:10.6f} | {a2*1000:10.6f} | {Z_in.real:12.2f} | {Z_in.imag:12.2f}")

        # Sweep frequency around 300 MHz to observe the parallel resonance
        # frequencies = np.linspace(290e6, 310e6, 11)
        
        # print("\nFreq (MHz) | Re(Z_in) [Ohms] | Im(Z_in) [Ohms]")
        # print("-" * 50)
        
        # for freq in frequencies:
        #     flat_z = compute_impedance_matrix(mesh, freq)
        #     Z = np.array(flat_z).reshape((N, N))
            
        #     I = la.solve(Z, V)
        #     Z_in = 1.0 / I[driven_dipole_idx]
            
        #     print(f"{freq/1e6:10.1f} | {Z_in.real:15.2f} | {Z_in.imag:15.2f}")
        # # end for
    # end for
    
    fig, ax = pyplot.subplots()
    
    for segment in mesh.segments:
        start_node = mesh.nodes[segment.start_idx]
        stop_node = mesh.nodes[segment.end_idx]
        
        ax.plot(
            [start_node.x, stop_node.x],
            [start_node.z, stop_node.z],
            marker='o',
        )
    # end for
    
    ax.grid()
    ax.legend()
    
    pyplot.show()

if __name__ == "__main__":
    main()