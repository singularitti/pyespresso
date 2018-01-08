#!/usr/bin/env python3
# created on Jan 6, 2018 at 08:31 by Qi Zhang

from typing import *

from json_tricks import load, loads


class JsonStr:
    """

    """

    def __init__(self, instr: Optional[str] = None, infile: Optional[str] = None):
        """


        :param instr: If this is given, *infile* argument will be ignored.
        :param infile: If *instream* is not given, this argument will be used.
        """
        if instr is None and infile is None:
            raise TypeError('instream and infile cannot be both None! You must specify one of them!')
        if isinstance(instr, str):
            self.instr: str = instr
            self.infile = None
        elif isinstance(infile, str):
            self.infile: str = infile
            self.instr = None
        else:
            raise TypeError('The type of one argument is wrong! They should all be str!')

    def content(self):
        if self.instr:
            return loads(self.instr)
        else:
            with open(self.infile, 'r') as f:
                return load(f)
