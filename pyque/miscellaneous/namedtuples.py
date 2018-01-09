#!/usr/bin/env python3
# created at Nov 28, 2017 2:08 PM by Qi Zhang
"""
:mod:`namedtuples` -- Toolkit for namedtuples
=============================================

.. module string
   :platform: Unix, Windows, Mac, Linux
   :synopsis: Analyze and convert namedtuple.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from collections import namedtuple


def is_named_tuple(obj: object) -> bool:
    """
    Check whether the object *obj* is a ``namedtuple``. This is referenced from
    `Pandas <https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L338-L362/>`_.

    :param obj: The object to check.
    :return: Whether the *obj* is a ``namedtuple``.

    .. doctest::

        >>> Point = namedtuple('Point', ['x', 'y'])
        >>> p = Point(1, 2)
        >>> is_named_tuple(p)
        True
        >>> is_named_tuple((1, 2))
        False
    """
    return isinstance(obj, tuple) and hasattr(obj, '_fields')
