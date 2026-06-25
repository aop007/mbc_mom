#!/usr/bin/env python

from typing import List

from . import geometry as geometry

def test_interface() -> None: ...

def compute_impedance_matrix(mesh: geometry.Mesh, frequency_hz: float) -> List[complex]: ...