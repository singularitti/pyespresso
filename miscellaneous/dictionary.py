#!/usr/bin/env python3
# created at Oct 22, 2017 12:05 AM by Qi Zhang


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts. Referenced from
    [here](https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression).

    :param dict_args: dictionaries
    :return: combined dictionary
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
