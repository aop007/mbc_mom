#!/usr/bin/env python



import numpy as np
import matplotlib.pyplot as plt
from mbc_mom.geometry import Mesh, Node, Segment
from mbc_mom import compute_impedance_matrix, compute_far_field, C

# 1. Initialize the Mesh
mesh = Mesh()

# 2. Define nodes along the Z-axis (length = 0.5m)
n0 = mesh.add_node(Node(0.0, 0.0, -0.25))
n1 = mesh.add_node(Node(0.0, 0.0,  0.00)) # Center node
n2 = mesh.add_node(Node(0.0, 0.0,  0.25))

# 3. Define the segments with a 1mm radius
radius = 0.001
mesh.add_segment(Segment(n0, n1, radius))
mesh.add_segment(Segment(n1, n2, radius))

# 4. Assemble the MBC basis functions
mesh.build_dipoles()

# 5. Validate the mesh geometry at our target frequency (300 MHz)
freq_hz = C  # Wavelength = 1 m
warnings = mesh.validate(freq_hz)
if warnings:
    print("Mesh validation warnings:", warnings)
    
# 6. Compute the Impedance Matrix [Z]
Z_list = compute_impedance_matrix(mesh, freq_hz)

# Reshape the flat list into a 2D NumPy array
# The dimension is N x N, where N is the number of dipoles
num_dipoles = len(mesh.dipoles)
Z_matrix = np.array(Z_list).reshape((num_dipoles, num_dipoles))

# 7. Define the excitation voltage vector [V]
# For a center-fed dipole, we apply 1V at the junction dipole
V = np.zeros(num_dipoles, dtype=complex)
V[0] = 1.0 + 0.0j # Assuming our single junction is at index 0

# 8. Solve for the currents [I] using standard linear algebra: [Z][I] = [V]
currents = np.linalg.solve(Z_matrix, V).tolist()

# 9. Define observation angles (Theta and Phi in radians)
# Elevation: 0 to 180 degrees
theta_array = np.linspace(0, np.pi, 181)
# Azimuth slice at Phi = 0
phi_array = np.array([np.pi / 2.0])                       

# 9.1 Build grids
thetas, phis = np.meshgrid(theta_array, phi_array, indexing='ij')
thetas = thetas.flatten()
phis = phis.flatten()

# 10. Compute the Far Field
E_far_field = compute_far_field(mesh, currents, freq_hz, list(zip(thetas, phis)))
E_far_field = np.array(E_far_field)
E_far_field = E_far_field.reshape(
    len(theta_array),
    len(phi_array),
)

print(f"Computed far field at {len(E_far_field)} angles.")

# 11. Display Normalized Elevation Pattern
norm_E_far_field = E_far_field / E_far_field.max()

fig, ax= plt.subplots(subplot_kw=dict(projection="polar"))

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

ax.plot(
    theta_array,
    norm_E_far_field[:, 0],
)

ax.set_rmin(0.0)
ax.set_rmax(1.0)
ax.set_rticks([0.25, 0.5, 0.75, 1.0])  # Fewer radial ticks
ax.grid(True)
ax.legend()

plt.show()