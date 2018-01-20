#!/usr/bin/env python3

import io
import sys
from typing import Optional, Iterator

from lazy_property import LazyProperty

# ========================================= What can be exported? =========================================
__all__ = ['TextStream']


class TextStream:
    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None, encoding: Optional[str] = None,
                 newline='\n'):
        """
        This is a general model for text streams. You can specify nothing in ``TextStream`` instance creation procedure,
        then data will be read from interactive input, you can press <Ctrl+d> to end this input.
        If you give it a string, with *newline* separator specified by ``'\n'`` or ``'\r\n'``,
        then the string will be parsed by this separator.
        If you give a file path in *infile*, and *instream* is not given or given but not a string, and if *infile* is
        a valid string for path, the *infile* will be read in ``stream_generator`` method.

        :param instream: If this is given, whatever *infile* you give will be ignored.
        :param infile: If *instream* is not given as a string, this argument will be used.
        :param encoding: Specifies the *infile*'s encoding when *instream* is not a string, and infile is given as a
            string. This keyword argument is just ``encoding`` argument for the builtin ``open`` function.
        :param newline: Specifies the *infile*'s newline character when *instream* is not a string,
            and infile is given as a string.
            This optional argument is just ``newline`` argument for the builtin ``open`` function.
            Or if a string or both first 2 arguments are ``None``, then this keyword argument is just the
            ``newline`` argument for ``io.StringIO``.
        """
        self.encoding = encoding
        self.newline = newline
        if instream is None and infile is None:
            self.instream: io.StringIO = io.StringIO(sys.stdin.read(), newline=newline)
            self.infile = infile  # ``None``
        elif isinstance(instream, str):
            # Here *infile* can be a string, ``None``, or any other type. Whatever it is, if *instream* is given,
            # we do not care.
            self.instream: io.StringIO = io.StringIO(instream, newline=newline)
            self.infile = None
        elif not isinstance(instream, str) and isinstance(infile, str):
            try:
                open(infile, 'r').close()
            except FileNotFoundError:
                raise FileNotFoundError("Input file '{0}' does not exist!".format(infile))
            # If *infile* is not a string, and *infile* is a string, then we care about *infile*.
            self.infile: str = infile
            self.instream = None
        else:
            raise TypeError('The type of at least one argument is wrong!')

    def stream_generator(self) -> Iterator[str]:
        """
        Create a generate that iterates the whole content of the file or string.

        :return: An iterator.
        """
        if self.instream:  # If *instream* is given and thus not ``None``.
            for line in self.instream:
                yield line
        else:  # If *infile* is given and thus not ``None``.
            with open(self.infile, 'r', encoding=self.encoding, newline=self.newline) as f:
                for line in f:
                    yield line

    @LazyProperty
    def contents(self) -> str:
        """
        Read the whole file or string, and return it.

        :return: The whole contents of the file or the string.
        """
        if isinstance(self.instream, io.StringIO):  # The first 2 possibilities in ``self.__init__``.
            return self.instream.getvalue()
        else:  # If *infile* is given but *instream* is not.
            with open(self.infile, 'r', encoding=self.encoding, newline=self.newline) as f:
                return f.read()
