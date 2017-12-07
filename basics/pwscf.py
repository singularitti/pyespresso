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
    def control_namelist(self):
        return self._control_namelist

    @control_namelist.setter
    def control_namelist(self, d: dict):
        self._control_namelist.update(d)

    @property
    def system_namelist(self):
        return self._system_namelist

    @system_namelist.setter
    def system_namelist(self, d: dict):
        self._system_namelist.update(d)

    @property
    def electrons_namelist(self):
        return self._electrons_namelist

    @electrons_namelist.setter
    def electrons_namelist(self, d: dict):
        self._electrons_namelist.update(d)

    @property
    def cell_parameters(self) -> np.ndarray:
        return self._cell_parameters

    @cell_parameters.setter
    def cell_parameters(self, cell_parameters: np.ndarray):
        self._cell_parameters = cell_parameters

    @property
    def atomic_species(self) -> List[NamedTuple]:
        return self._atomic_species

    @atomic_species.setter
    def atomic_species(self, atomic_species: List[NamedTuple]):
        self._atomic_species = atomic_species

    @property
    def atomic_positions(self) -> List[NamedTuple]:
        return self._atomic_positions

    @atomic_positions.setter
    def atomic_positions(self, atomic_positions: List[NamedTuple]):
        self._atomic_positions = atomic_positions

    @property
    def k_points(self) -> NamedTuple:
        return self._k_points

    @k_points.setter
    def k_points(self, k_points: NamedTuple):
        self._k_points = k_points

    @property
    def atomic_positions_option(self) -> str:
        if self._atomic_positions_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        return self._atomic_positions_option

    @atomic_positions_option.setter
    def atomic_positions_option(self, option: str):
        if self._atomic_positions_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        self._atomic_positions_option = option

    @property
    def k_points_option(self) -> str:
        return self._k_points_option

    @k_points_option.setter
    def k_points_option(self, option: str):
        self._k_points_option = option

    def write_to_file(self, out_file: str):
        k_points = self._k_points
        with open(out_file, 'w') as f:
            f.write("&CONTROL\n")
            for k, v in self.control_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&SYSTEM\n")
            for k, v in self.system_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&ELECTRONS\n")
            for k, v in self.electrons_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\nCELL_PARAMETERS\n")
            f.write(
                re.sub("[\[\]]", ' ',
                       np.array2string(self._cell_parameters,
                                       formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
            )
            f.write("\nATOMIC_SPECIES\n")
            for row in self._atomic_species:
                f.write(' '.join(row))
                f.write("\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self._atomic_positions_option))
            for row in self._atomic_positions:
                f.write(' '.join(row))
                f.write("\n")
            f.write("K_POINTS {{ {0} }}\n".format(self._k_points_option))
            f.write(' '.join(map(str, (k_points.grid + k_points.offsets))))
        print("Object is written to file '{0}'!".format(out_file))
