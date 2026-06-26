#!/usr/bin/env python

from typing import List

from . import geometry as geometry

def test_interface() -> None: ...

def compute_impedance_matrix(mesh: geometry.Mesh, frequency_hz: float) -> List[complex]: ...

def compute_far_field(mesh: geometry.Mesh, currents: List[complex], freq_hz: float, thetas: List[float], phis: List[float]) -> List[float]: ...
