#!/usr/bin/env python3
# created at Dec 3, 2017 2:35 AM by Qi Zhang

from numbers import Number

import numpy as np


def is_number(obj: object):
    """
    Check if the object is a number.
    Referenced from [Pandas](https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L27-L48).

    Examples
    --------
    >>> is_number(np.float64(1))
    True
    >>> is_number("foo")
    False

    :param obj: The object to check.
    :return: Whether `obj` is a number or not.
    """
    return isinstance(obj, (Number, np.number))
