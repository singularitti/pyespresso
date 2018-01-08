#!/usr/bin/env python3
# created at Nov 28, 2017 2:08 PM by Qi Zhang

from collections import namedtuple


def is_named_tuple(obj: object) -> bool:
    """
    Check if the object is a named tuple. This is from
    [Pandas library](https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L338-L362).

    Examples
    --------
    >>> Point = namedtuple('Point', ['x', 'y'])
    >>> p = Point(1, 2)
    >>> is_named_tuple(p)
    True
    >>> is_named_tuple((1, 2))
    False

    :param obj: The object to check.
    :return: Whether `obj` is a named tuple.
    """
    return isinstance(obj, tuple) and hasattr(obj, '_fields')
