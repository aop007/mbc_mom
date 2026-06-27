#!/usr/bin/env python

from typing import List, Optional

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
    def length(self, nodes: List['Node']) -> float: ...

class Dipole:
    seg1_idx: int
    seg2_idx: int
    junction_idx: int
    mbc_offset: float
    is_monopole: bool
    def __init__(self, seg1_idx: int, seg2_idx: int, junction_idx: int, mbc_offset: float) -> None: ...
    
class GroundPlane:
    is_pec: bool
    sigma: float
    eps_r: float
    use_sommerfeld: bool

class Mesh:
    nodes: List[Node]
    segments: List[Segment]
    dipoles: List[Dipole]
    ground_plane: Optional[GroundPlane]
    
    def __init__(self) -> None: ...
    def add_node(self, node: Node) -> int: ...
    def add_segment(self, segment: Segment) -> None: ...
    def build_dipoles(self) -> None: ...
    def set_pec_ground(self) -> None: ...
    def set_real_ground(self, sigma: float, eps_r: float, use_sommerfeld: bool) -> None: ...
    def validate(self, freq_hz: float) -> List[str]: ...