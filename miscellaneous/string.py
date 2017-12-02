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

    Examples
    --------
    >>> strs_to_ints(['1', '1.0', '-0.2'])
    [1, 1, 0]

    :param strs: a list of string
    :return: a list of converted integers
    """
    return _strs_to_(strs, lambda x: int(float(x)))


def strs_to_floats(strs: List[str]) -> List[float]:
    """
    Convert a list of strings to a list of floats.

    Examples
    --------
    >>> strs_to_floats(['1', '1.0', '-0.2'])
    [1.0, 1.0, -0.2]

    :param strs: a list of string
    :return: a list of converted floats
    """
    return _strs_to_(strs, float)


def str_to_double_precision_float(s: str) -> float:
    """
    Double precision float in Fortran file will have form 'x.ydz' or 'x.yDz', this cannot be convert directly to float
    by Python `float` function, so I wrote this function to help conversion. For example,

    Examples
    --------
    >>> str_to_double_precision_float('1d-82')
    1e-82
    >>> str_to_double_precision_float('1.0D-82')
    1e-82
    >>> str_to_double_precision_float('0.8D234')
    8e+233
    >>> str_to_double_precision_float('.8d234')
    8e+233

    :param s: a string denoting a double precision number
    :return: a Python floating point number
    """
    first, second, exponential = re.match("(-?\d*)\.?(-?\d*)d(-?\d+)", s, re.IGNORECASE).groups()
    return float(first + '.' + second + 'e' + exponential)


def str_to_float(s: str) -> float:
    """
    Convert a string to corresponding single or double precision scientific number.

    Examples
    --------
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
            return str_to_double_precision_float(s)
        except AttributeError:
            raise ValueError("The string '{0}' does not corresponds to a double precision number!".format(s))
    else:
        return float(s)


def match_one_string(pattern: str, s: str, *args):
    """
    Make sure you know only none or one string will be matched! If you are not sure, use `match_one_pattern` instead.

    Examples
    --------
    >>> p = "\d+"
    >>> s = "abc 123 def"
    >>> match_one_string(p, s, int)
    123
    >>> print(match_one_string(p, "abc"))
    Pattern "\d+" not found, or more than one found in string abc!
    None
    >>> print(match_one_string(p, "abc 123 def 456"))
    Pattern "\d+" not found, or more than one found in string abc 123 def 456!
    None

    :param pattern:
    :param s:
    :param args:
    :return:
    """
    try:
        match, = re.findall(pattern, s)  # `match` is either an empty list or a list of string.
        if len(args) == 0:  # If no wrapper argument is given, return directly the matched string
            return match
        elif len(args) == 1:  # If wrapper argument is given, i.e., not empty, then apply wrapper to the match
            wrapper, = args
            return wrapper(match)
        else:
            raise TypeError('Multiple wrappers are given! Only one should be given!')
    except ValueError:
        print("Pattern \"{0}\" not found, or more than one found in string {1}!".format(pattern, s))


def match_one_pattern(pattern: str, s: str, *args: Optional[Callable], **flags):
    """
    Find a pattern in a certain string. If found and a wrapper is given, then return the wrapped matched-string; if no
    wrapper is given, return the pure matched string. If no match is found, return None.

    Examples
    --------
    >>> p = "\d+"
    >>> s = "abc 123 def 456"
    >>> match_one_pattern(p, s)
    ['123', '456']
    >>> match_one_pattern(p, s, int)
    [123, 456]
    >>> match_one_pattern(p, "abc 123 def")
    ['123']
    >>> print(match_one_pattern('s', 'abc'))
    Pattern "s" not found in string abc!
    None
    >>> match_one_pattern('s', 'Ssa', flags=re.IGNORECASE)
    ['S', 's']

    :param pattern: a pattern, can be a string or a regular expression
    :param s: a string
    :param args: at most 1 argument can be given
    :param flags: the same flags as `re.findall`'s
    :return:
    """
    match: Optional[List[str]] = re.findall(pattern, s,
                                            **flags)  # `match` is either an empty list or a list of strings.
    if match:
        if len(args) == 0:  # If no wrapper argument is given, return directly the matched string
            return match
        elif len(args) == 1:  # If wrapper argument is given, i.e., not empty, then apply wrapper to the match
            wrapper, = args
            return [wrapper(m) for m in match]
        else:
            raise TypeError('Multiple wrappers are given! Only one should be given!')
    else:  # If no match is found
        print("Pattern \"{0}\" not found in string {1}!".format(pattern, s))
        return None


def is_string_like(obj: object) -> bool:
    """
    Check if the object is a string.

    Examples
    --------
    >>> is_string_like("foo")
    True
    >>> is_string_like(1)
    False

    :param obj: The object to check.
    :return: Whether `obj` is a string or not.
    """
    return isinstance(obj, str)
