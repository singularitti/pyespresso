#!/usr/bin/env python3
# created at Dec 3, 2017 10:46 PM by Qi Zhang

import os
import re

import numpy as np
from lazy_property import LazyWritableProperty


class PHononStandaradInput:
    def __init__(self):
        self.__name__ = 'PHononStandaradInput'
        self._INPUTPH_namelist: dict = {}
        self._single_q_point = None
        self._q_points = None

    @LazyWritableProperty
    def title(self) -> str:
        pass

    @property
    def INPUTPH_namelist(self):
        return self._INPUTPH_namelist

    @INPUTPH_namelist.setter
    def INPUTPH_namelist(self, d: dict):
        self._INPUTPH_namelist.update(d)

    @LazyWritableProperty
    def single_q_point(self) -> np.ndarray:
        pass

    @LazyWritableProperty
    def q_points(self) -> np.ndarray:
        pass

    def write_to_file(self, out_file: str):
        with open(out_file, 'w') as f:
            f.write(self.title)
            f.write("/\n&INPUTPH\n")
            for k, v in self._INPUTPH_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            if self._single_q_point:
                f.write(' '.join(self._single_q_point))
                f.write("\n")
            elif self._q_points:
                f.write(
                    re.sub("[\[\]]", ' ',
                           np.array2string(self._q_points,
                                           formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
                )
            else:
                print("Object '{0}' is written to file {1}!".format(self.__name__, os.path.abspath(out_file)))
