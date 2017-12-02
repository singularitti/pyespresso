#!/usr/bin/env python3
# created at Dec 1, 2017 2:25 PM by Qi Zhang

import os


def is_file_like(obj) -> bool:
    """
    Check if the object is a file-like object.
    For objects to be considered file-like, they must be an iterator and
    have either a `read` and/or `write` method as an attribute.
    Note: file-like objects must be iterable, but iterable objects need not be file-like.

    Examples
    --------
    >>> print(os.path.abspath(__file__))
    >>> is_file_like(os.path.abspath(__file__))
    True
    >>> is_file_like([1, 2, 3])
    False

    :param obj: The object to check.
    :return: Whether `obj` has file-like properties.
    """
    if not (hasattr(obj, 'read') or hasattr(obj, 'write')):
        return False

    if not hasattr(obj, "__iter__"):
        return False

    return True
