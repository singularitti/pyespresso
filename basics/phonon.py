#!/usr/bin/env python3
# created at Dec 3, 2017 10:46 PM by Qi Zhang

import numpy as np
from basics.lazy import CachedProperty


class PHononStandaradInput:
    def __init__(self):
        self._title: str = ''
        self._INPUTPH_namelist: dict = {}
        self._phonon_wavevector: np.ndarray = np.empty(3)
        self._q_points: np.ndarray = np.empty(4)

    @CachedProperty
    def title(self):
        return self._title

    @title.setter
    def title(self, title: str):
        self._title = title

    @CachedProperty
    def INPUTPH_namelist(self):
        return self._INPUTPH_namelist

    @INPUTPH_namelist.setter
    def INPUTPH_namelist(self, d: dict):
        self._INPUTPH_namelist.update(d)

    @CachedProperty
    def phonon_wavevector(self):
        return self._phonon_wavevector

    @phonon_wavevector.setter
    def phonon_wavevector(self, wavevector: np.ndarray):
        self._phonon_wavevector = wavevector

    @CachedProperty
    def q_points(self):
        return self._q_points

    @q_points.setter
    def q_points(self, q_points: np.ndarray):
        self._q_points = q_points

    def write_to_file(self, out_file: str):
        with open(out_file, 'w') as f:
            f.write(self._title)
            f.write("/\n&INPUTPH\n")
            for k, v in self._INPUTPH_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
