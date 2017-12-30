#!/usr/bin/env python3
# created at Dec 3, 2017 10:46 PM by Qi Zhang

import os
import re

import numpy as np


class PHononStandaradInput:
    def __init__(self):
        self.__name__ = 'PHononStandaradInput'
        self.__title__: str = ''
        self._INPUTPH_namelist: dict = {}
        self._single_q_point = None
        self._q_points = None

    @property
    def title(self):
        return self.__title__

    @title.setter
    def title(self, title: str):
        self.__title__ = title

    @property
    def INPUTPH_namelist(self):
        return self._INPUTPH_namelist

    @INPUTPH_namelist.setter
    def INPUTPH_namelist(self, d: dict):
        self._INPUTPH_namelist.update(d)

    @property
    def single_q_point(self):
        return self._single_q_point

    @single_q_point.setter
    def single_q_point(self, q_point: np.ndarray):
        self._single_q_point = q_point

    @property
    def q_points(self):
        return self._q_points

    @q_points.setter
    def q_points(self, q_points: np.ndarray):
        self._q_points = q_points

    def write_to_file(self, out_file: str):
        with open(out_file, 'w') as f:
            f.write(self.__title__)
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
