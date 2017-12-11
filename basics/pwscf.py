#!/usr/bin/env python3
# created at Nov 21, 2017 2:13 AM by Qi Zhang

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
        raise TypeError('{0} is not a pw input object!'.format(obj))


class PWStandardInput:
    def __init__(self):
        self._control_namelist: Dict[str, Any] = {}
        self._system_namelist: Dict[str, Any] = {}
        self._electrons_namelist: Dict[str, Any] = {}
        self._cell_parameters: np.ndarray = np.empty((3, 3))
        self._k_points: NamedTuple = None
        self._atomic_species: List[NamedTuple] = []
        self._atomic_positions: List[NamedTuple] = []
        self._cell_parameters_option: str = ''
        self._atomic_positions_option: str = ''
        self._k_points_option: str = ''

    @property
    def CONTROL(self):
        return self._control_namelist

    @CONTROL.setter
    def CONTROL(self, d: dict):
        self._control_namelist.update(d)

    @property
    def SYSTEM(self):
        return self._system_namelist

    @SYSTEM.setter
    def SYSTEM(self, d: dict):
        self._system_namelist.update(d)

    @property
    def ELECTRONS(self):
        return self._electrons_namelist

    @ELECTRONS.setter
    def ELECTRONS(self, d: dict):
        self._electrons_namelist.update(d)

    @property
    def CELL_PARAMETERS(self) -> np.ndarray:
        return self._cell_parameters

    @CELL_PARAMETERS.setter
    def CELL_PARAMETERS(self, cell_parameters: np.ndarray):
        self._cell_parameters = cell_parameters

    @property
    def ATOMIC_SPECIES(self) -> List[NamedTuple]:
        return self._atomic_species

    @ATOMIC_SPECIES.setter
    def ATOMIC_SPECIES(self, atomic_species: List[NamedTuple]):
        self._atomic_species = atomic_species

    @property
    def ATOMIC_POSITIONS(self) -> List[NamedTuple]:
        return self._atomic_positions

    @ATOMIC_POSITIONS.setter
    def ATOMIC_POSITIONS(self, atomic_positions: List[NamedTuple]):
        self._atomic_positions = atomic_positions

    @property
    def K_POINTS(self) -> NamedTuple:
        return self._k_points

    @K_POINTS.setter
    def K_POINTS(self, k_points: NamedTuple):
        self._k_points = k_points

    @property
    def ATOMIC_POSITIONS_option(self) -> str:
        if self._atomic_positions_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        return self._atomic_positions_option

    @ATOMIC_POSITIONS_option.setter
    def ATOMIC_POSITIONS_option(self, option: str):
        if self._atomic_positions_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        self._atomic_positions_option = option

    @property
    def K_POINTS_option(self) -> str:
        return self._k_points_option

    @K_POINTS_option.setter
    def K_POINTS_option(self, option: str):
        self._k_points_option = option

    def write_to_file(self, out_file: str):
        k_points = self._k_points
        with open(out_file, 'w') as f:
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
                       np.array2string(self._cell_parameters,
                                       formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
            )
            f.write("\nATOMIC_SPECIES\n")
            for row in self._atomic_species:
                f.write(' '.join(map(str, row)))
                f.write("\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self._atomic_positions_option))
            for row in self._atomic_positions:
                f.write(' '.join(row))
                f.write("\n")
            f.write("K_POINTS {{ {0} }}\n".format(self._k_points_option))
            f.write(' '.join(map(str, (k_points.grid + k_points.offsets))))
        print("Object is written to file '{0}'!".format(out_file))
