#!/usr/bin/env python3
# created at Dec 28, 2017 2:45 PM by Qi Zhang
"""
:mod:`sets` -- Toolkit for sets
===============================

.. module sets
   :platform: Unix, Windows, Mac, Linux
   :synopsis: Analyze and convert sets.
.. moduleauthor Qi Zhang <qz2280@columbia.edu>
"""


def add_elements_to_set(s: set, *args) -> set:
    """
    Add one or more elements in a set *s*.

    :param s: The set to be added elements.
    :param args: The element(s) to be added.
    :return: The modified set *s*.

    .. doctest::

        >>> p = set()
        >>> add_elements_to_set(p, 1)
        {1}
        >>> q = set()
        >>> add_elements_to_set(q)
        set()
        >>> s = set()
        >>> add_elements_to_set(s, 1, 2, 3)
        {1, 2, 3}
    """
    s.update(set(*args))
    return s


def remove_elements_from_set(s: set, *args) -> set:
    """
    Remove one or more elements in a set *s*.

    :param s: The set to be removed elements.
    :param args: The element(s) to be removed.
    :return: The modified set *s*.

    .. doctest::

        >>> s = {1, 2, 3, 4, 5}
        >>> remove_elements_from_set(s, 1)
        {2, 3, 4, 5}
        >>> remove_elements_from_set(s, 2, 3, 4)
        {5}
        >>> remove_elements_from_set(s)
        {5}
    """
    for _ in args:
        s.remove(_)
    return s
