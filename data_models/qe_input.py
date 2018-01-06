#!/usr/bin/env python3
# created at Nov 21, 2017 2:13 AM by Qi Zhang

import os
import re
import warnings
from collections import namedtuple
from typing import *

import numpy as np
from lazy_property import LazyWritableProperty

from meta.descriptors import LabeledDescriptor, MetaDescriptorOwner
from miscellaneous.path_generator import path_generator

# ========================================= type alias =========================================
KPoints = NamedTuple('KPoints', [('grid', int), ('offsets', int)])
AtomicSpecies = NamedTuple('AtomicSpecies', [('name', str), ('mass', float), ('pseudopotential', str)])
AtomicPosition = NamedTuple('AtomicPosition', [('name', str), ('position', np.ndarray)])

# ========================================= define useful data structures =========================================
KPoints: KPoints = namedtuple('KPoints', ['grid', 'offsets'])

AtomicSpecies: AtomicSpecies = namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])
AtomicSpecies.__doc__ = \
    """Note that the word 'species' serves as singular and plural both.
    So here though suffixed with a 's', it is for one atom and thus is singular."""

AtomicPosition: AtomicPosition = namedtuple('AtomicPosition', ['name', 'position'])


# ========================================= define useful functions =========================================
def is_pw_input(obj: object):
    if all(hasattr(obj, attr) for attr in ['CONTROL', 'SYSTEM', 'ELECTRONS', 'CELL_PARAMETERS']):
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
        raise TypeError("{0} is not a {1} object!".format(obj, type(obj).__name__))


# ========================================= define useful data structures =========================================
class _Option(LabeledDescriptor):
    def __set__(self, instance, new_option: str):
        if new_option == 'alat':
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        else:
            instance.__dict__[self.label] = new_option


# ========================================= most important data structures =========================================
class PWscfStandardInput(metaclass=MetaDescriptorOwner):
    atomic_positions_option = _Option()

    def __init__(self):
        self._CONTROL: Dict[str, Any] = {}
        self._SYSTEM: Dict[str, Any] = {}
        self._ELECTRON: Dict[str, Any] = {}

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

    @LazyWritableProperty
    def cell_parameters(self) -> np.ndarray:
        pass

    @LazyWritableProperty
    def atomic_species(self) -> List[NamedTuple]:
        pass

    @LazyWritableProperty
    def atomic_positions(self) -> List[NamedTuple]:
        pass

    @LazyWritableProperty
    def k_points(self) -> NamedTuple:
        pass

    @LazyWritableProperty
    def k_points_option(self) -> str:
        pass

    def to_text_file(self, outfile: str, path_prefix: Optional[str] = '') -> None:
        outfile_path = path_generator(outfile, path_prefix)

        k_points: KPoints = self.k_points

        with open(outfile_path, 'w') as f:
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
            f.write(re.sub("[\[\]]", ' ', np.array2string(self.cell_parameters,
                                                          formatter={'float_kind': lambda x: "{:20.10f}".format(x)})))
            f.write("\nATOMIC_SPECIES\n")
            for row in self.atomic_species:
                f.write(' '.join(map(str, row)) + "\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self.atomic_positions_option))
            for row in self.atomic_positions:
                f.write(' '.join(row) + "\n")
            f.write("K_POINTS {{ {0} }}\n".format(self.k_points_option))
            f.write(' '.join(map(str, (k_points.grid + k_points.offsets))))

        print("Object '{0}' is written to file {1}!".format(self.__name__, os.path.abspath(outfile_path)))


class PHononStandardInput:
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

    def write_to_file(self, outfile: str, path_prefix: Optional[str] = '') -> None:
        outfile_path = path_generator(outfile, path_prefix)

        with open(outfile_path, 'w') as f:
            f.write(self.title)
            f.write("/\n&INPUTPH\n")
            for k, v in self._INPUTPH_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            if self._single_q_point:
                f.write(' '.join(self._single_q_point))
                f.write("\n")
            elif self._q_points:
                f.write(re.sub("[\[\]]", ' ', np.array2string(self._q_points,
                                                              formatter={
                                                                  'float_kind': lambda x: "{:20.10f}".format(x)})))
            else:
                print("Object '{0}' is written to file {1}!".format(self.__name__, os.path.abspath(outfile_path)))
