#!/usr/bin/env python3
# created at Oct 22, 2017 12:05 AM by Qi Zhang
"""
:mod:`dicts` -- Toolkit for dictionaries
========================================

.. module dicts
   :platform: Unix, Windows, Mac, Linux
   :synopsis: Analyze and convert dictionaries.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""


def merge_dicts(*dicts):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts. Referenced from
    `here <https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression/>`_.
    Note this will have the values merged, i.e., the latter dict's values overwrite those from the former.

    .. note::

        Existing keys will be overwritten by the last ``dict``!

    :param dicts: dictionaries
    :return: a flattened, combined dictionary

    .. doctest::

        >>> a = {'a': 0, 'b': 1}
        >>> b = {'c': 2, 'd': 3}
        >>> c = {'e': 4, 'f': 5}
        >>> d = {'a': 1, 'f': 2}
        >>> merge_dicts(a, b, c, d)
        {'a': 1, 'b': 1, 'c': 2, 'd': 3, 'e': 4, 'f': 2}
    """
    result = {}
    for dictionary in dicts:
        result.update(dictionary)  # Note ``update`` overwrites existing keys!
    return result


def is_dict_like(obj: object) -> bool:
    """
    Check if the object is dict-like. Referenced from
    `Pandas <https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L365-L394/>`_.

    :param obj: The object to check.
    :return: Whether *obj* has dict-like properties.

    .. doctest::

        >>> is_dict_like({1: 2})
        True
        >>> is_dict_like([1, 2, 3])
        False
    """
    return hasattr(obj, '__getitem__') and hasattr(obj, 'keys')


def list_public_attributes(obj: object):
    return [k for k, v in vars(obj).items() if not (k.startswith('_') or callable(v))]
