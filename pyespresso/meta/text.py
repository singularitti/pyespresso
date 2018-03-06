#!/usr/bin/env python3
"""
:mod:`text` -- Text Data Model
==============================

.. module text
   :platform: Unix, Windows, Mac, Linux
   :synopsis:
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import pathlib
from io import StringIO
from sys import stdin
from typing import Optional, Iterator, Tuple, Union

from lazy_property import LazyProperty

from pyespresso.tools.os_path import is_path_exists_or_creatable

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

    :param inp: Input, can be a string, an ``StringIO`` object, or ``None`` (which means read from standard input).
        If the *inp* is a valid path for the system, the file the *inp* directs will be used.

    .. [#f] Referenced from `here <https://docs.microsoft.com/en-us/cpp/c-runtime-library/text-and-binary-streams>`_.
    """

    def __init__(self, inp: Union[str, StringIO, None] = None, **kwargs):
        self.__infile_path = None

        if inp is None:
            self.__stream = StringIO(stdin.read())
        elif isinstance(inp, str):
            if pathlib.Path(inp).expanduser().is_file():
                self.__infile_path = inp
                with open(inp, **kwargs) as f:
                    self.__stream = StringIO(f.read())
            else:
                self.__stream: StringIO = StringIO(inp)
        elif isinstance(inp, StringIO):
            self.__stream = inp
        else:
            raise TypeError("The type '{0}' of *inp* argument is not supported!".format(type(inp)))

    @LazyProperty
    def stream(self) -> StringIO:
        """
        Read-only property.

        :return:
        """
        return self.__stream

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
        stream = self.__stream
        stream.seek(0)
        for line in stream:
            yield line

    def generator_telling_position(self) -> Iterator[Tuple[str, int]]:
        """
        Create a generate that iterates the whole content of the file or string, and also tells which offset is now.

        :return: An iterator.
        """
        stream = self.__stream
        stream.seek(0)
        for line in stream:
            yield line, stream.tell()

    def generator_starts_from(self, offset, whence: Optional[int] = 0) -> Iterator[str]:
        """

        :param offset:
        :param whence:
        :return:
        """
        stream = self.__stream
        stream.seek(offset, whence)
        for line in stream:
            yield line

    @LazyProperty
    def content(self) -> str:
        """
        Read the whole file or string, and return it.

        :return: The whole contents of the file or the string.
        """
        return self.__stream.getvalue()

    def to_file(self, filename: str) -> None:
        """
        Write this object to file named *filename*.

        :param filename:
        :return: ``None``
        """
        import shutil

        if not is_path_exists_or_creatable(filename):
            raise OSError('PATH is not reachable or exists!')

        with open(filename, 'w') as f:
            self.__stream.seek(0)
            shutil.copyfileobj(self.__stream, f)
