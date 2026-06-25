#!/usr/bin/env python

import mbc_mom

def main():
    # Initialize the mesh container
    mesh = mbc_mom.Mesh()

    # Create a step-radius dipole along the Z-axis
    n0 = mesh.add_node(mbc_mom.Node(0.0, 0.0, -0.5))
    n1 = mesh.add_node(mbc_mom.Node(0.0, 0.0, 0.0))  # The Junction Node
    n2 = mesh.add_node(mbc_mom.Node(0.0, 0.0, 0.5))

    # Add segments with different radii meeting at n1
    seg1 = mbc_mom.Segment(n0, n1, radius=0.001)
    seg2 = mbc_mom.Segment(n1, n2, radius=0.005)
    
    mesh.add_segment(seg1)
    mesh.add_segment(seg2)

    print(f"Mesh created with {len(mesh.nodes)} nodes and {len(mesh.segments)} segments.")
    print(f"Segment 1 length: {seg1.length(mesh.nodes):.3f} m")
# end main()


if __name__ == "__main__":
    main()
# end if