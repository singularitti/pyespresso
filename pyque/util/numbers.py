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
