#!/usr/bin/env python3
"""
:mod:`numbers` -- Toolkit for numbers
=====================================

.. module numbers
   :platform: Unix, Windows, Mac, Linux
   :synopsis: Analyze and convert namedtuple.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from numbers import Number
from typing import *

import numpy as np


def is_number(obj: object):
    """
    Check if the object is a number. Referenced from
    `Pandas <https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L27-L48/>`_.

    :param obj: The object to check.
    :return: Whether *obj* is a number or not.

    .. doctest::

        >>> is_number(float(5))
        True
        >>> is_number(np.float64(1))
        True
        >>> is_number("foo")
        False
    """
    return isinstance(obj, (Number, np.number))


def is_integer(obj: object) -> bool:
    """
    Check if the object is an integer.

    :param obj: The object to check.
    :return: Whether *obj* is an integer or not.

    .. doctest::

        >>> is_integer(2.0)
        True
        >>> is_integer(2)
        True
        >>> is_integer('s')
        False
    """
    if is_number(obj):
        if isinstance(obj, int):
            return True
        elif isinstance(obj, float):
            return obj.is_integer()
        else:
            return False
    else:
        import warnings
        warnings.warn("Only numbers can be tested if they are integers!", stacklevel=2)
        return False


def is_float(obj: object) -> bool:
    if is_number(obj):
        if isinstance(obj, float):
            return True
        else:
            return False
    else:
        return False


def all_integer_like(iterable: Iterable[object]):
    return all(is_integer(_) for _ in iterable)


def all_float_like(iterable: Iterable[object]):
    return all(is_float(_) for _ in iterable)
