import numpy as np
import scipy.linalg as la
from typing import List, Tuple, Optional
from mbc_mom import (
    compute_impedance_matrix, 
    get_incident_eval_points, 
    compute_incident_v_matrix
)
from mbc_mom.geometry import Mesh, Node, Segment

# --- 1. MOCKING YOUR EXTERNAL SYSTEM ---

class MockCoordinate:
    def __init__(self, lon: float, lat: float, elevation: float):
        self.lon = lon
        self.lat = lat
        self.elevation = elevation

class MockRadiator:
    def __init__(self, coord: MockCoordinate):
        self.coord = coord

class MockPropagationModel:
    """Mocks your PropagationModel returning a plane wave arriving from the Y-axis."""
    def get_field_strength(self, tx: MockRadiator, xs: np.ndarray, ys: np.ndarray, zs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        # Assume a 1 V/m plane wave polarized along the Z-axis, traveling along -Y
        # E(y) = E0 * exp(j * k * y)
        freq_hz = 300e6
        k = (2.0 * np.pi * freq_hz) / 299792458.0
        
        # Calculate Phase shifts across the grid
        phase = k * ys
        E_mag = 1.0
        
        Ex = np.zeros_like(xs, dtype=np.complex128)
        Ey = np.zeros_like(ys, dtype=np.complex128)
        Ez = (E_mag * np.exp(1j * phase)).astype(np.complex128)
        
        # We don't need H for the MoM [V] matrix, so we return dummy zeros
        Hx, Hy, Hz = np.zeros_like(Ex), np.zeros_like(Ey), np.zeros_like(Ez)
        
        return (Ex, Ey, Ez), (Hx, Hy, Hz)

# --- 2. THE PIPELINE IMPLEMENTATION ---

def build_receiving_dipole(mesh: Mesh, num_segments: int):
    """Builds a vertical half-wave dipole along the Z-axis."""
    z_coords = np.linspace(-0.25, 0.25, num_segments + 1)
    nodes = [mesh.add_node(Node(0.0, 0.0, float(z))) for z in z_coords]
    for i in range(num_segments):
        mesh.add_segment(Segment(nodes[i], nodes[i+1], 0.001))
    return nodes[num_segments // 2]

def main():
    print("--- MBC Solver: External Incident Field Scattering ---")
    mesh = Mesh()
    center_node = build_receiving_dipole(mesh, 21)
    mesh.build_dipoles()
    
    N = len(mesh.dipoles)
    freq_hz = 300e6
    points_per_seg = 7
    
    # Setup your external model environment
    analysis_coord = MockCoordinate(-71.0, 46.0, 100.0)
    distant_tx = MockRadiator(MockCoordinate(-71.1, 46.1, 150.0))
    prop_model = MockPropagationModel()
    
    # 1. Evaluate [Z] Matrix (Self-impedance of the receiving antenna)
    print("1. Computing MoM [Z] Matrix...")
    Z = np.array(compute_impedance_matrix(mesh, freq_hz)).reshape((N, N))
    
    # 2. Extract local coordinates for the incident field
    print(f"2. Extracting Integration Points ({points_per_seg} per segment)...")
    xs_flat, ys_flat, zs_flat = get_incident_eval_points(mesh, points_per_seg)
    
    # 3. Query the Propagation Model
    # Here is where you would map (xs_flat, ys_flat, zs_flat) to Geodetic.
    # For this script, we just pass the raw local ENU coordinates directly.
    print("3. Querying PropagationModel for E-fields...")
    (Ex, Ey, Ez), _ = prop_model.get_field_strength(
        distant_tx, 
        np.array(xs_flat), 
        np.array(ys_flat), 
        np.array(zs_flat)
    )
    
    # 4. Integrate Incident Field into MoM [V] Matrix
    print("4. Integrating [V] Matrix in Rust...")
    V = np.array(compute_incident_v_matrix(
        mesh, freq_hz, 
        Ex.tolist(), Ey.tolist(), Ez.tolist(), 
        points_per_seg
    ))
    
    # 5. Solve for Induced Currents
    print("5. Solving I = Z^-1 V...")
    I = la.solve(Z, V)
    
    # Analysis
    center_idx = next(i for i, d in enumerate(mesh.dipoles) if d.junction_idx == center_node)
    received_current = I[center_idx]
    
    print("\n--- Results ---")
    print(f"Incident Field Magnitude: 1.0 V/m")
    print(f"Induced Current at Center: {abs(received_current) * 1000:.2f} mA")
    print(f"Current Phase: {np.angle(received_current, deg=True):.2f} degrees")

if __name__ == "__main__":
    main()