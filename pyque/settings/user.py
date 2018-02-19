#!/usr/bin/env python3
"""
:mod:`mod` -- title
========================================

.. module mod
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import collections
import pathlib

import json_tricks

# ========================================= What can be exported? =========================================
__all__ = ['Settings', 'from_json']


class Settings(collections.ChainMap):
    def to_json_file(self, filename: str):
        if not filename.endswith('.json'):
            filename += '.json'
        with open(filename, 'w') as f:
            json_tricks.dump(self.maps, f)

    def to_json_string(self):
        json_tricks.dumps(self.maps)


def from_json(inp: str) -> Settings:
    if pathlib.Path(inp).is_file():
        return Settings(json_tricks.load(inp))
    else:
        return Settings(json_tricks.loads(inp))
