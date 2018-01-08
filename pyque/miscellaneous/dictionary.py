#!/usr/bin/env python3
# created at Oct 22, 2017 12:05 AM by Qi Zhang


def merge_dicts(*dicts):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts. Referenced from
    [here](https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression).
    Note this will have the values merged, i.e., the latter dict's values overwrite those from the former.

    Examples
    --------
    >>> a = {'a': 0, 'b': 1}
    >>> b = {'c': 2, 'd': 3}
    >>> c = {'e': 4, 'f': 5}
    >>> d = {'a': 1, 'f': 2}
    >>> merge_dicts(a, b, c, d)
    {'a': 1, 'b': 1, 'c': 2, 'd': 3, 'e': 4, 'f': 2}

    :param dicts: dictionaries
    :return: a flattened, combined dictionary
    """
    result = {}
    for dictionary in dicts:
        result.update(dictionary)  # `update` overwrites existing keys
    return result


def is_dict_like(obj: object) -> bool:
    """
    Check if the object is dict-like. Referenced from
    [Pandas](https://github.com/pandas-dev/pandas/blob/v0.21.0/pandas/core/dtypes/inference.py#L365-L394).

    Examples
    --------
    >>> is_dict_like({1: 2})
    True
    >>> is_dict_like([1, 2, 3])
    False

    :param obj: The object to check.
    :return: Whether `obj` has dict-like properties.
    """
    return hasattr(obj, '__getitem__') and hasattr(obj, 'keys')
