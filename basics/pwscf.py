#!/usr/bin/env python3
# created at Nov 21, 2017 2:13 AM by Qi Zhang

import os
import re
import warnings
from collections import namedtuple
from typing import *

import numpy as np

KPoints = namedtuple('KPoints', ['grid', 'offsets'])
AtomicSpecies = namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])
AtomicPositions = namedtuple('AtomicPositions', ['name', 'positions'])


def is_pw_input(obj: object):
    if all(hasattr(obj, attr) for attr in
           ['CONTROL', 'SYSTEM', 'ELECTRONS', 'CELL_PARAMETERS']):
        return True
    else:
        return False


def print_pw_input(obj: object):
    if is_pw_input(obj):
        try:
            from beeprint import pp
            pp(obj)
        except ModuleNotFoundError:
            print(obj)
    else:
        raise TypeError("{0} is not a 'PWStandardInput' object!".format(obj))


class PWStandardInput:
    def __init__(self):
        self.__name__ = 'PWStandardInput'
        self._CONTROL: Dict[str, Any] = {}
        self._SYSTEM: Dict[str, Any] = {}
        self._ELECTRON: Dict[str, Any] = {}
        self._CELL_PARAMETERS: np.ndarray = np.empty((3, 3))
        self._K_POINTS: NamedTuple = None
        self._ATOMIC_SPECIES: List[NamedTuple] = []
        self._ATOMIC_POSITIONS: List[NamedTuple] = []
        self._CELL_PARAMETERS_option: str = ''
        self._ATOMIC_POSITIONS_option: str = ''
        self._K_POINTS_option: str = ''

    @property
    def CONTROL(self):
        return self._CONTROL

    @CONTROL.setter
    def CONTROL(self, d: dict):
        self._CONTROL.update(d)

    @property
    def SYSTEM(self):
        return self._SYSTEM

    @SYSTEM.setter
    def SYSTEM(self, d: dict):
        self._SYSTEM.update(d)

    @property
    def ELECTRONS(self):
        return self._ELECTRON

    @ELECTRONS.setter
    def ELECTRONS(self, d: dict):
        self._ELECTRON.update(d)

    @property
    def CELL_PARAMETERS(self) -> np.ndarray:
        return self._CELL_PARAMETERS

    @CELL_PARAMETERS.setter
    def CELL_PARAMETERS(self, cell_parameters: np.ndarray):
        self._CELL_PARAMETERS = cell_parameters

    @property
    def ATOMIC_SPECIES(self) -> List[NamedTuple]:
        return self._ATOMIC_SPECIES

    @ATOMIC_SPECIES.setter
    def ATOMIC_SPECIES(self, atomic_species: List[NamedTuple]):
        self._ATOMIC_SPECIES = atomic_species

    @property
    def ATOMIC_POSITIONS(self) -> List[NamedTuple]:
        return self._ATOMIC_POSITIONS

    @ATOMIC_POSITIONS.setter
    def ATOMIC_POSITIONS(self, atomic_positions: List[NamedTuple]):
        self._ATOMIC_POSITIONS = atomic_positions

    @property
    def K_POINTS(self) -> NamedTuple:
        return self._K_POINTS

    @K_POINTS.setter
    def K_POINTS(self, k_points: KPoints):
        self._K_POINTS = k_points

    @property
    def ATOMIC_POSITIONS_option(self) -> str:
        if self._ATOMIC_POSITIONS_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        return self._ATOMIC_POSITIONS_option

    @ATOMIC_POSITIONS_option.setter
    def ATOMIC_POSITIONS_option(self, option: str):
        if self._ATOMIC_POSITIONS_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        self._ATOMIC_POSITIONS_option = option

    @property
    def K_POINTS_option(self) -> str:
        return self._K_POINTS_option

    @K_POINTS_option.setter
    def K_POINTS_option(self, option: str):
        self._K_POINTS_option = option

    def write_to_file(self, output_path: str):
        k_points: KPoints = self._K_POINTS
        with open(output_path, 'w') as f:
            f.write("&CONTROL\n")
            for k, v in self.CONTROL.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&SYSTEM\n")
            for k, v in self.SYSTEM.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&ELECTRONS\n")
            for k, v in self.ELECTRONS.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\nCELL_PARAMETERS\n")
            f.write(
                re.sub("[\[\]]", ' ',
                       np.array2string(self._CELL_PARAMETERS,
                                       formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
            )
            f.write("\nATOMIC_SPECIES\n")
            for row in self._ATOMIC_SPECIES:
                f.write(' '.join(map(str, row)))
                f.write("\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self._ATOMIC_POSITIONS_option))
            for row in self._ATOMIC_POSITIONS:
                f.write(' '.join(row))
                f.write("\n")
            f.write("K_POINTS {{ {0} }}\n".format(self._K_POINTS_option))
            f.write(' '.join(map(str, (k_points.grid + k_points.offsets))))
        print("Object '{0}' is written to file {1}!".format(self.__name__, os.path.abspath(output_path)))
