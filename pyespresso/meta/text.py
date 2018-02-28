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
from typing import Optional, Iterator, Tuple, Union

from lazy_property import LazyProperty

# ========================================= What can be exported? =========================================
__all__ = ['TextStream']


class TextStream:
    """
    This is a general model for text streams.
    A text stream consists of one or more lines of text that can be written to a text-oriented display
    so that they can be read. When reading from a text stream, the program reads a *newline* at
    the end of each line. When writing to a text stream, the program writes a *newline* to signal the
    end of a line. [#f]_

    You can specify nothing in ``TextStream`` instance creation procedure,
    then data will be read from interactive input, you can press <Ctrl+D> to end this input.
    If you give it a string, with *newline* separator specified by ``'\\n'`` or ``'\\r\\n'``,
    then the string will be parsed by this separator.
    If you give a file path in *inp*, and if it is valid, then the file will be read.

    :param inp:
    :param newline: This optional argument is just the ``newline`` argument for the builtin ``open`` function
        or the argument for ``StringIO``.
    :param kwargs:

    .. [#f] Referenced from `here <https://docs.microsoft.com/en-us/cpp/c-runtime-library/text-and-binary-streams>`_.
    """

    def __init__(self, inp: Union[str, io.StringIO, None] = None, newline: Optional[str] = None, **kwargs):
        self.newline = newline
        self.__open_kwargs: dict = kwargs

        self.__infile_path = None

        if inp is None:
            self.__stream: io.StringIO = io.StringIO(sys.stdin.read(), newline=newline)
        elif isinstance(inp, str):
            if pathlib.Path(inp).expanduser().is_file():
                self.__infile_path = inp
                self.__stream = None
            else:
                self.__stream: io.StringIO = io.StringIO(inp, newline=newline)
        elif isinstance(inp, io.StringIO):
            self.__stream = inp
        else:
            raise TypeError("The type '{0}' of *inp* argument is not supported!".format(type(inp)))

    @LazyProperty
    def stream(self) -> io.StringIO:
        """
        Read-only property.

        :return:
        """
        if self.__infile_path is None:
            return self.__stream
        else:
            with open(self.__infile_path, newline=self.newline, **self.__open_kwargs) as f:
                return io.StringIO(f.read())

    @LazyProperty
    def infile_path(self) -> Optional[pathlib.PurePath]:
        """
        Read-only property.

        :return:
        """
        if not self.__infile_path:
            return None
        else:
            return pathlib.Path(self.__infile_path).expanduser()

    def generator(self) -> Iterator[str]:
        """
        Create a generate that iterates the whole content of the file or string.

        :return: An iterator.
        """
        if self.__infile_path is None:
            stream = self.__stream
            for line in stream:
                yield line
        else:
            with open(self.__infile_path, newline=self.newline, **self.__open_kwargs) as f:
                for line in f:
                    yield line

    def generator_telling_position(self) -> Iterator[Tuple[str, int]]:
        """
        Create a generate that iterates the whole content of the file or string, and also tells which offset is now.

        :return: An iterator.
        """
        if self.__infile_path is None:
            stream = self.__stream
            for line in stream:
                yield line, stream.tell()
        else:
            with open(self.__infile_path, newline=self.newline, **self.__open_kwargs) as f:
                for line in iter(f.readline, ''):
                    yield line, f.tell()

    def generator_starts_from(self, offset, whence: Optional[int] = 0) -> Iterator[str]:
        """


        :param offset:
        :param whence:
        :return:
        """
        if self.__infile_path is None:
            stream = self.__stream
            stream.seek(offset, whence)
            for line in stream:
                yield line
        else:
            with open(self.__infile_path, newline=self.newline, **self.__open_kwargs) as f:
                f.seek(offset, whence)
                for line in f:
                    yield line

    @LazyProperty
    def content(self) -> str:
        """
        Read the whole file or string, and return it.

        :return: The whole contents of the file or the string.
        """
        if self.__infile_path is None:
            return self.__stream.getvalue()
        else:
            with open(self.__infile_path, newline=self.newline, **self.__open_kwargs) as f:
                return f.read()

    def to_file(self, filename: str) -> None:
        """
        Write this object to file named *filename*.

        :param filename:
        :return: ``None``
        """
        if not pathlib.Path(filename).expanduser().exists():
            raise FileNotFoundError()

        with open(filename, 'w') as f:
            self.__stream.seek(0)
            shutil.copyfileobj(self.__stream, f)
