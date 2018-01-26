#!/usr/bin/env python3
"""
:mod:`mod` --
========================================

.. module mod
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import re
from typing import *

import numpy as np

from pyque.meta.text import TextStream
from pyque.tools.strings import strings_to_integers, strings_to_floats


class FreqLexer:
    def __init__(self, inp: Optional[str]):
        self.inp = inp
        self.text_stream = TextStream(inp)
        self.commenters = '!'
        self.nbnd = None
        self.nks = None

    def lex_output(self) -> Tuple[List[float], ...]:
        generator: Iterator[str] = self.text_stream.stream_generator()
        for line in generator:
            if '&plot' in line:
                match = re.search("nbnd=\s*(\d+).*nks=\s*(\d+)", line, flags=re.IGNORECASE)
                if match is None:
                    raise RuntimeError('Something went wrong!')
                else:
                    try:
                        self.nbnd, self.nks = strings_to_integers(match.groups())
                    except ValueError:
                        raise
                    break
        bands = []
        qs = []
        for line in generator:
            line_stripped = line.strip()
            if not line_stripped or line_stripped.startswith(self.commenters):  # empty line or comment line
                continue
            _: List[float] = strings_to_floats(line_stripped.split())
            if len(_) == self.nbnd:
                bands.append(_)
            elif len(_) == 3:
                qs.append(_)
            else:
                raise ValueError('Unknown value read!')
        return bands, qs


class FreqParser:
    def __init__(self, lexer: FreqLexer):
        if isinstance(lexer, FreqLexer):
            self.lexer = lexer
        else:
            raise TypeError("The argument you gave for *lexer* is not an instance of 'FreqLexer'!")
        self.nbnd = None
        self.nks = None
        self.bands = None
        self.qs = None

    def check(self) -> bool:
        nbnd = self.lexer.nbnd
        nks = self.lexer.nks
        try:
            bands, qs = self.lexer.lex_output()
        except ValueError:
            raise
        self.bands = np.array(bands)
        self.qs = np.array(qs)
        if self.bands.shape == (self.lexer.nks, self.lexer.nbnd) and self.qs.shape == (self.lexer.nks, 3):
            self.nbnd = self.lexer.nbnd
            self.nks = self.lexer.nks
            return True
        else:
            if self.bands.shape != (nks, nbnd):
                print('The shape of bands is incorrect!')
            else:
                print('The shape of q-points is incorrect!')
            return False

    def separate_bands(self) -> Optional[List[np.ndarray]]:
        if not self.check():
            return None
        else:
            bands = []
            for i in range(self.nbnd):
                bands.append(self.bands[:, i])
            return bands


def auto_parse(inp: Optional[str]) -> Tuple[Optional[List[np.ndarray]], np.ndarray]:
    parser = FreqParser(FreqLexer(inp))
    return parser.separate_bands(), parser.qs
