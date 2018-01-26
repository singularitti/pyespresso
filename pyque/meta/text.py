#!/usr/bin/env python3
"""
:mod:`text` -- Text Data Model
========================================

.. module text
   :platform: Unix, Windows, Mac, Linux
   :synopsis:
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import io
import pathlib
import shutil
import sys
from typing import Optional, Iterator, Tuple

# ========================================= What can be exported? =========================================
__all__ = ['TextStream']


class TextStream:
    def __init__(self, inp: Optional[str] = None, encoding: Optional[str] = None, newline: Optional[str] = None):
        """
        This is a general model for text streams. You can specify nothing in ``TextStream`` instance creation procedure,
        then data will be read from interactive input, you can press <Ctrl+d> to end this input.
        If you give it a string, with *newline* separator specified by ``'\n'`` or ``'\r\n'``,
        then the string will be parsed by this separator.
        If you give a file path in *infile*, and *instream* is not given or given but not a string, and if *infile* is
        a valid string for path, the *infile* will be read in ``stream_generator`` method.

        :param inp: If this is given, whatever *infile* you give will be ignored.
        :param encoding: Specifies the *infile*'s encoding when *instream* is not a string, and infile is given as a
            string. This keyword argument is just ``encoding`` argument for the builtin ``open`` function.
        :param newline: This optional argument is just the ``newline`` argument for the builtin ``open`` function
            or an argument for ``StringIO``.
        """
        self.encoding = encoding
        self.newline = newline
        self.infile = None
        if inp is None:
            self.instream: io.StringIO = io.StringIO(sys.stdin.read(), newline=newline)
        elif isinstance(inp, str):
            if pathlib.Path(inp).is_file():
                self.infile = inp
            else:
                self.instream: io.StringIO = io.StringIO(inp, newline=newline)
        else:
            raise TypeError("The type of *inp* argument '{0}' is wrong! It should be a string!".format(type(inp)))

    def stream_generator(self) -> Iterator[str]:
        """
        Create a generate that iterates the whole content of the file or string.

        :return: An iterator.
        """
        if self.infile is None:  # The default *self.infile* is set to ``None``.
            instream = self.instream
            for line in instream:
                yield line
        else:
            with open(self.infile, 'r', encoding=self.encoding, newline=self.newline) as f:
                for line in f:
                    yield line

    def stream_generator_with_position(self) -> Iterator[Tuple[str, int]]:
        if self.infile is None:  # The default *self.infile* is set to ``None``.
            instream = self.instream
            for line in instream:
                yield line, instream.tell()
        else:
            with open(self.infile, 'r', encoding=self.encoding, newline=self.newline) as f:
                for line in iter(f.readline, ''):
                    yield line, f.tell()

    @property
    def contents(self) -> str:
        """
        Read the whole file or string, and return it.

        :return: The whole contents of the file or the string.
        """
        if self.infile is None:
            return self.instream.getvalue()
        else:  # If *infile* is given but *instream* is not.
            with open(self.infile, 'r', encoding=self.encoding, newline=self.newline) as f:
                return f.read()

    def to_text_file(self, file: str):
        with open(file, 'w') as f:
            self.instream.seek(0)
            shutil.copyfileobj(self.instream, f)
