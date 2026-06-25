from typing import List

class Node:
    x: float
    y: float
    z: float
    def __init__(self, x: float, y: float, z: float) -> None: ...

class Segment:
    start_idx: int
    end_idx: int
    radius: float
    def __init__(self, start_idx: int, end_idx: int, radius: float) -> None: ...
    def length(self, nodes: List[Node]) -> float: ...

class Dipole:
    seg1_idx: int
    seg2_idx: int
    junction_idx: int
    mbc_offset: float
    def __init__(self, seg1_idx: int, seg2_idx: int, junction_idx: int, mbc_offset: float) -> None: ...

class Mesh:
    nodes: List[Node]
    segments: List[Segment]
    dipoles: List[Dipole]
    def __init__(self) -> None: ...
    def add_node(self, node: Node) -> int: ...
    def add_segment(self, segment: Segment) -> None: ...
    def build_dipoles(self) -> None: ...