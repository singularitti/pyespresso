#!/usr/bin/env python3
# Created on Jan 5, 2018 at 20:02 by Qi Zhang

from io import StringIO
from typing import Optional, Iterator

from lazy_property import LazyProperty


class TextStream:
    def __init__(self, instr: Optional[str] = None, infile: Optional[str] = None):
        """


        :param instr: If this is given, `infile` argument will be ignored.
        :param infile: If `instream` is not given, this argument will be used.
        """
        if instr is None and infile is None:
            raise TypeError('instream and infile cannot be both None! You must specify one of them!')
        if isinstance(instr, str):
            self.instream: str = instr
            self.infile = None
        elif isinstance(infile, str):
            self.infile: str = infile
            self.instream = None
        else:
            raise TypeError('The type of one argument is wrong! They should all be str!')

    def stream_generator(self) -> Iterator[str]:
        if self.instream:  # If `instream` is given and thus not `None`
            for line in self.instream.split("\n"):
                yield line
        else:  # If `infile` is given and thus not `None`
            with open(self.infile, 'r') as f:
                for line in f:
                    yield line

    @LazyProperty
    def to_string_io(self):
        if self.instream:
            return StringIO(self.instream)
        else:
            with open(self.infile, 'r') as f:
                return StringIO(f.read())
