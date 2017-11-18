#!/usr/bin/env python3
# created at Oct 22, 2017 12:05 AM by Qi Zhang


def merge_dicts(*dicts):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts. Referenced from
    [here](https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression).
    Note this will have the values merged, i.e., the latter dict's values overwrite those from the former.

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
