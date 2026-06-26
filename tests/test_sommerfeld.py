import numpy as np
import scipy.linalg as la
from mbc_mom import compute_impedance_matrix
from mbc_mom.geometry import Mesh, Node, Segment

def build_horizontal_dipole(mesh: Mesh, num_segments: int, height: float):
    """Builds a 0.5m horizontal dipole along the X-axis at a specific Z height."""
    x_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = [mesh.add_node(Node(float(x), 0.0, height)) for x in x_coords]
    
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.001))
        
    return nodes[num_segments // 2]

def solve_antenna(use_sommerfeld: bool) -> complex:
    mesh = Mesh()
    
    # Wet Earth: Moderately high permittivity and conductivity
    eps_r = 15.0
    sigma = 0.01
    mesh.set_real_ground(sigma, eps_r, use_sommerfeld)
    
    # Place the dipole just 5 cm above the dirt
    driven_node_idx = build_horizontal_dipole(mesh, 20, height=0.05)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = 300e6
    
    # Solve MoM
    Z = np.array(compute_impedance_matrix(mesh, freq_hz)).reshape((N, N))
    V = np.zeros(N, dtype=np.complex128)
    
    driven_dipole_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == driven_node_idx)
    V[driven_dipole_idx] = 1.0 + 0.0j

    I = la.solve(Z, V)
    
    return 1.0 / I[driven_dipole_idx]

def main():
    print("--- MBC Solver: Ground Wave Interaction Benchmark ---")
    print("Antenna: 300 MHz Horizontal Dipole")
    print("Height:  5 cm above Wet Earth (eps_r=15, sigma=0.01)\n")
    
    print("1. Solving with Near-Field RCA (Approximation)...")
    z_rca = solve_antenna(use_sommerfeld=False)
    print(f"   Z_in = {z_rca.real:.2f} + {z_rca.imag:.2f}j Ohms\n")
    
    print("2. Generating LUT and Solving with Rigorous Sommerfeld Integrals...")
    z_somm = solve_antenna(use_sommerfeld=True)
    print(f"   Z_in = {z_somm.real:.2f} + {z_somm.imag:.2f}j Ohms\n")
    
    print("--- Analysis ---")
    diff_real = abs(z_rca.real - z_somm.real)
    diff_imag = abs(z_rca.imag - z_somm.imag)
    print(f"Real Impedance Error in RCA: {diff_real:.2f} Ohms")
    print(f"Imag Reactance Error in RCA: {diff_imag:.2f} Ohms")

if __name__ == "__main__":
    main()