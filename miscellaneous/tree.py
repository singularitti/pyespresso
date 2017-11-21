#!/usr/bin/env python3
# created at Nov 21, 2017 4:31 PM by Qi Zhang
"""
Referenced from [here](https://stackoverflow.com/a/43237270/3260253).
"""

from typing import *


class Tree(dict):
    """
    A tree implementation using python's autovivification feature.
    """

    def __init__(self, data: Optional[Any] = {}):
        """
        Cast a (nested) dict to a (nested) `Tree` class

        :param data: can be anything, normal type, dict, and even `Tree`
        """
        super().__init__()

        for k, v in data.items():
            if isinstance(v, dict):
                # A recursive setup of a new `Tree`,
                # `type(self)(v)` is new a `Tree` with `v` as value and `k` as key, thus a nested `Tree` is constructed
                self[k] = type(self)(v)
            else:
                self[k] = data

    def __missing__(self, k):
        # `type(self)()` creates a new empty `Tree`, and use it as the value of the current tree's key `k`, and return
        # this value.
        v = self[k] = type(self)()
        return v


if __name__ == '__main__':
    t = Tree()
    t['a']['b']['x'] = 't.a.b.x'
    t['a']['b']['y'] = 't.a.b.y'
    t['a']['c']['x'] = 't.a.c.x'
    t['a']['d'] = 't.a.d'
    t['b'] = 't.b'
    print('t is:\n{0}'.format(t))
