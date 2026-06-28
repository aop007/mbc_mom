#!/usr/bin/env python

from typing import Optional

import click
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

from mbc_mom import compute_impedance_matrix, C
from mbc_mom.geometry import Mesh, Node, Segment


def build_balun(
    freq_Hz: float,
    guy_radius_m: float,
    anchor_angle_deg: float,
    balun_anchor_dist_m: float,
    balun_length_m: float,
    balun_radius_m: float,
    balun_branchs: int,
    balun_source_dist_m: float,
    segment_resolution_m: float,
) -> Mesh:
    mesh = Mesh()
    
    raise NotImplementedError()
# end 