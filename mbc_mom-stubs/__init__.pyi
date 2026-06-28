#!/usr/bin/env python

from typing import List, Tuple

from . import geometry as geometry

C: float
EPS_0: float
J: float
ETA: float

def test_interface() -> None: ...

def compute_impedance_matrix(mesh: geometry.Mesh, frequency_hz: float) -> List[complex]: ...

def compute_far_field(mesh: geometry.Mesh, currents: List[complex], freq_hz: float, thetas_phis: List[tuple[float, float]]) -> List[float]: ...

def compute_near_field(
    mesh: geometry.Mesh, 
    currents: List[complex], 
    freq_hz: float, 
    xs: List[float], 
    ys: List[float], 
    zs: List[float]
) -> Tuple[List[complex], List[complex], List[complex], List[complex], List[complex], List[complex]]: ...

def get_incident_eval_points(
    mesh: geometry.Mesh, 
    points_per_seg: int = 7
) -> Tuple[List[float], List[float], List[float]]: ...

def compute_incident_v_matrix(
    mesh: geometry.Mesh,
    freq_hz: float,
    ex: List[complex],
    ey: List[complex],
    ez: List[complex],
    points_per_seg: int = 7
) -> List[complex]: ...