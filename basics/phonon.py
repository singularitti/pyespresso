#!/usr/bin/env python3
# created at Dec 3, 2017 10:46 PM by Qi Zhang

import numpy as np


class PHononStandaradInput:
    def __init__(self):
        self._INPUTPH_namelist: dict = {}
        self._phonon_wavevector: np.ndarray = np.empty(3)
        self._qPointsSpecs: np.ndarray = np.empty(4)

    @property
    def INPUTPH_namelist(self):
        return self._INPUTPH_namelist

    @INPUTPH_namelist.setter
    def INPUTPH_namelist(self, d: dict):
        self._INPUTPH_namelist.update(d)

    @property
    def phonon_wavevector(self):
        return self._phonon_wavevector

    @phonon_wavevector.setter
    def phonon_wavevector(self, wavevector: np.ndarray):
        self._phonon_wavevector = wavevector
