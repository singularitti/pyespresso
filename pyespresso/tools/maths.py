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


def derivative1d(a: Union[np.ndarray, List[float]], b: Union[np.ndarray, List[float]]):
    """
    Given 2 array-like objects, compute the first-order derivative for *a* W.R.T *b*.

    :param a: An array-like object, same length as *b*.
    :param b: An array-like object, same length as *a*.
    :return: A numpy array which has the same length as *a* and *b*.

    .. doctest::

        >>> a = np.array([1, 2, 3, 5, 6])
        >>> b = np.array([0.5, 1.5, 1, 2, 3])
        >>> derivative1d(a, b)
        array([1. , 4. , 6. , 1.5, 1. ])
    """
    return np.gradient(a) / np.gradient(b)
