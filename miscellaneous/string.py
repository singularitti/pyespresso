#!/usr/bin/env python3
# created at Nov 17, 2017 3:56 PM by Qi Zhang

import re
from typing import *


def _strs_to_(strs: List[str], wrapper: Callable) -> List:
    """
    Convert a list of strings to a list of certain form, specified by `wrapper`.

    :param strs: a list of string
    :param wrapper: a function that converts your string
    :return: type undefined, but specified by `to_type`
    """
    return list(map(wrapper, strs))


def strs_to_ints(strs: List[str]) -> List[int]:
    """
    Convert a list of strings to a list of integers.

    >>> strs_to_ints(['1', '1.0', '-0.2'])
    [1, 1, 0]

    :param strs: a list of string
    :return: a list of converted integers
    """
    return _strs_to_(strs, lambda x: int(float(x)))


def strs_to_floats(strs: List[str]) -> List[float]:
    """
    Convert a list of strings to a list of floats.

    >>> strs_to_floats(['1', '1.0', '-0.2'])
    [1.0, 1.0, -0.2]

    :param strs: a list of string
    :return: a list of converted floats
    """
    return _strs_to_(strs, float)


def str_to_double_scientific_float(s: str) -> float:
    """
    Double precision float in Fortran file will have form 'x.ydz' or 'x.yDz', this cannot be convert directly to float
    by Python `float` function, so I wrote this function to help conversion. For example,

    >>> str_to_double_scientific_float('1d-82')
    1e-82
    >>> str_to_double_scientific_float('1.0D-82')
    1e-82
    >>> str_to_double_scientific_float('0.8D234')
    8e+233
    >>> str_to_double_scientific_float('.8d234')
    8e+233

    :param s: a string denoting a double precision number
    :return: a Python floating point number
    """
    first, second, exponential = re.match("(-?\d*)\.?(-?\d*)d(-?\d+)", s, re.IGNORECASE).groups()
    return float(first + '.' + second + 'e' + exponential)


def str_to_float(s: str) -> float:
    """
    Convert a string to corresponding single or double precision scientific number.

    >>> str_to_float('1.0D-5')
    1e-05
    >>> str_to_float('1Dx')
    Traceback (most recent call last):
        ...
    ValueError: The string '1Dx' does not corresponds to a double precision number!
    >>> str_to_float('.8d234')
    8e+233
    >>> str_to_float('0.1')
    0.1

    :param s: a string could be '0.1', '1e-5', '1.0D-5', or any other validated number
    :return: a float or raise an error
    """
    if 'D' in s.upper():  # Possible double precision number
        try:
            return str_to_double_scientific_float(s)
        except AttributeError:
            raise ValueError("The string '{0}' does not corresponds to a double precision number!".format(s))
    else:
        return float(s)
