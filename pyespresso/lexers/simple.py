#!/usr/bin/env python3
"""
:mod:`mod` -- title
========================================

.. module mod
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import io
import re
from typing import *

from pyespresso.meta.text import TextStream

# ========================================= What can be exported? =========================================
__all__ = ['SimpleLexer']


# ================================= These are some type aliases or type definitions. =================================

# ========================================= define useful data structures =========================================
class SimpleLexer:
    def __init__(self, inp: Union[str, io.StringIO, None] = None, **kwargs):
        self.text_stream = TextStream(inp, **kwargs)

    def match_one_string(self, pattern: str, flags: Optional[int] = 0) -> Optional[str]:
        """
        This method matches a pattern which exists only once in the file.

        :param pattern: a regular expression that you want to match
        :param flags:
        :return: Determined by the `wrapper`, the value you want to grep out from pyespresso.the file.
        """
        text: str = self.text_stream.content
        regex: Pattern = re.compile(pattern, flags)
        match: Optional[Match] = regex.search(text)

        if not match:  # If no match is found
            print("Pattern {0} is not found!".format(pattern))
            return None

        return match.group(0)

    def _match_to_float(self, pattern: str) -> Optional[float]:
        _: Optional[str] = self.match_one_string(pattern)

        if not _:  # If no match is found
            return None

        return float(_)

    def _match_to_int(self, pattern: str) -> Optional[int]:
        _: Optional[str] = self.match_one_string(pattern)

        if not _:  # If no match is found
            return None

        return int(_)

    def _read_n_columns(self, n: int) -> List[List[str]]:
        """
        Read any number of columns (if exists) separated by spaces.

        :param n: an integer specifies the number of columns you want to read
        :return: a list of `n` columns that contain the contents of each column
        """
        n_columns = [[] for _ in range(n)]
        for line in self.text_stream.content:
            if not line.split():  # If line is ''
                continue
            else:
                sp = line.split()
                try:
                    for i, column in enumerate(n_columns):
                        column.append(sp[i])
                except IndexError:
                    print('You want more columns than you have! Check your file at line {0}!'.format(line))
        return n_columns

    def read_two_columns(self) -> List[List[str]]:
        """
        This method reads 2 columns from pyespresso.a file. The are specified by spaces.

        :return: two columns of the file
        """
        return self._read_n_columns(2)

    def _read_one_column_as_keys(self, col_index: int, func: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param col_index: the index of the column that you want to make it as keys.
        :param func: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        keys = []
        values = []
        # Add utf-8 support because we may use special characters.
        for line in self.text_stream.generator():
            sp = line.split()
            keys.append(sp[col_index])
            del sp[col_index]  # Remove the indexing column
            values.append(func(sp))
        return dict(zip(keys, values))
