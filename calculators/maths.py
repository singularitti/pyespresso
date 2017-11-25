#!/usr/bin/env python3
# created at Aug 26, 2017 9:04 PM by Qi Zhang

from typing import *

import numpy as np


def compute_3d_distance(point1: List[float], point2: List[float]) -> np.float64:
    """
    This function computes 3D Euclidean distance between 2 points.

    >>> np.testing.assert_equal(compute_3d_distance([0, 0, 0], [1, 1, 1]), np.sqrt(3))

    :param point1: point 1
    :param point2: point 2
    :return: scalar distance
    """
    point1 = np.array(point1)
    point2 = np.array(point2)
    return np.sqrt(((point2 - point1) ** 2).sum())
