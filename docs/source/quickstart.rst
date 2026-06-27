Quickstart: Half-Wave Dipole
============================

This guide demonstrates how to build and simulate a simple center-fed half-wave dipole using ``mbc_mom``. We will construct the mesh, compute the impedance matrix, solve for the antenna currents, and evaluate the far-field pattern.

Building the Mesh
-----------------

First, we define the geometry. A half-wave dipole at 300 MHz has a wavelength of roughly 1 meter, so we will create a dipole that is 0.5 meters long, oriented along the Z-axis.

To resolve the center feed-point, we define three nodes and two segments.

.. code-block:: python

    import numpy as np
    from mbc_mom.geometry import Mesh, Node, Segment
    from mbc_mom import compute_impedance_matrix, compute_far_field

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
    freq_hz = 300e6
    warnings = mesh.validate(freq_hz)
    if warnings:
        print("Mesh validation warnings:", warnings)


Solving the System
------------------

With the geometry defined, we pass the mesh to the highly parallelized Rust backend to compute the MoM impedance matrix. 

.. code-block:: python

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


Evaluating the Far Field
------------------------

Finally, we use the solved currents to compute the far-field radiation pattern.

.. code-block:: python

    # 9. Define observation angles (Theta and Phi in radians)
    thetas = np.linspace(0, np.pi, 181).tolist() # Elevation: 0 to 180 degrees
    phis = [0.0] * 181                           # Azimuth slice at Phi = 0

    # 10. Compute the Far Field
    E_far_field = compute_far_field(mesh, currents, freq_hz, thetas, phis)

    print(f"Computed far field at {len(E_far_field)} angles.")