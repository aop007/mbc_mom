#!/usr/bin/env python


import numpy as np
from mbc_mom.geometry import Mesh, Node, Segment
from mbc_mom import compute_impedance_matrix

mesh = Mesh()

# Step-radius junction
n0 = mesh.add_node(Node(0.0, 0.0, -0.5))
n1 = mesh.add_node(Node(0.0, 0.0, 0.0))  # Junction node
n2 = mesh.add_node(Node(0.0, 0.0, 0.5))

mesh.add_segment(Segment(n0, n1, radius=0.001))
mesh.add_segment(Segment(n1, n2, radius=0.005))

# Add parallel VDipole
n3 = mesh.add_node(Node(1.0, 0.0, -0.5))
n4 = mesh.add_node(Node(1.0, 0.0, 0.0))  # Junction node
n5 = mesh.add_node(Node(1.0, 0.0, 0.5))

mesh.add_segment(Segment(n3, n4, radius=0.001))
mesh.add_segment(Segment(n4, n5, radius=0.005))

mesh.build_dipoles()

# Get the flat matrix from Rust
flat_z = compute_impedance_matrix(mesh, 100e6) # 100 MHz

# Reshape into a 2D NxN NumPy array
N = len(mesh.dipoles)
Z_matrix = np.array(flat_z).reshape((N, N))

print(f"Matrix shape: {Z_matrix.shape}")
print("Exact Reciprocity Check (Z_01 == Z_10):", Z_matrix[0,1] == Z_matrix[1,0])