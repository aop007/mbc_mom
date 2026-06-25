#!/usr/bin/env python

from mbc_mom import Mesh, Node, Segment

def test_mbc_topology():
    mesh = Mesh()

    # Step-radius junction
    n0 = mesh.add_node(Node(0.0, 0.0, -0.5))
    n1 = mesh.add_node(Node(0.0, 0.0, 0.0))  # Junction node
    n2 = mesh.add_node(Node(0.0, 0.0, 0.5))

    mesh.add_segment(Segment(n0, n1, radius=0.001))
    mesh.add_segment(Segment(n1, n2, radius=0.005))

    # Generate the interaction dipoles
    mesh.build_dipoles()

    print(f"Total Segments: {len(mesh.segments)}")
    print(f"Total Dipoles: {len(mesh.dipoles)}")
    for i, d in enumerate(mesh.dipoles):
        print(f"Dipole {i}: Segments ({d.seg1_idx}, {d.seg2_idx}) at Junction Node {d.junction_idx}")
        print(f"  -> MBC Offset Applied: {d.mbc_offset:.4f} m")
    # end for
# end test_mbc_topology()

if __name__ == "__main__":
    test_mbc_topology()
# end if