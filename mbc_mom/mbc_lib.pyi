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

class Mesh:
    nodes: List[Node]
    segments: List[Segment]
    def __init__(self) -> None: ...
    def add_node(self, node: Node) -> int: ...
    def add_segment(self, segment: Segment) -> None: ...