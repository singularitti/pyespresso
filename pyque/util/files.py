#!/usr/bin/env python3
"""
:mod:`files` -- Dealing with files
==================================

.. module sets
   :platform: Unix, Windows, Mac, Linux
   :synopsis:
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""


def is_file_like(obj: object) -> bool:
    """
    Check if the object is a file-like object.
    For objects to be considered file-like, they must be an iterator and
    have either a ``read`` and/or ``write`` method as an attribute.
    The code is referenced from `pandas.api.types.is_file_like
    <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.api.types.is_file_like.html#pandas.api.types.is_file_like/>`_.

    .. note::

        File-like objects must be iterable, but iterable objects need not be file-like.

    :param obj: The object to check.
    :return: Whether *obj* has file-like attributes.

    .. testsetup::

        from io import StringIO

    .. doctest::

        >>> from io import StringIO
        >>> is_file_like(StringIO('data'))
        True
        >>> is_file_like([1, 2, 3])
        False
    """
    if not (hasattr(obj, 'read') or hasattr(obj, 'write')):
        return False

    if not hasattr(obj, "__iter__"):
        return False

    return True
