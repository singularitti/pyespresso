#!/usr/bin/env python3

import os
import re
import warnings
from collections import namedtuple
from typing import *

import numpy as np
from json_tricks import dump, dumps
from lazy_property import LazyWritableProperty

from pyque.meta.card import LazyCard
from pyque.meta.descriptors import LabeledDescriptor, DescriptorOwnerMeta
from pyque.meta.namelist import LazyNamelist, NamelistDict
from pyque.util.path_generators import path_generator

# ========================================= What can be exported? =========================================
__all__ = ['is_pw_input', 'print_pw_input', 'PWscfInput', 'SCFInput', 'VCRelaxInput',
           'PHononStandardInput', 'AtomicSpecies', 'AtomicPosition', 'KPoints']

# ========================================= type alias =========================================
KPoints = NamedTuple('KPoints', [('grid', int), ('offsets', int)])
AtomicSpecies = NamedTuple('AtomicSpecies', [('name', str), ('mass', float), ('pseudopotential', str)])
AtomicPosition = NamedTuple('AtomicPosition', [('name', str), ('position', np.ndarray)])

# ========================================= define useful data structures =========================================
KPoints: KPoints = namedtuple('KPoints', ['grid', 'offsets'])

AtomicSpecies: AtomicSpecies = namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])
AtomicSpecies.__doc__ = """\
Note that the word 'species' serves as singular and plural both. 
So here though suffixed with a 's', it is for one atom and thus is singular."""

AtomicPosition: AtomicPosition = namedtuple('AtomicPosition', ['name', 'position'])


# ========================================= variables declaration =========================================


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
class _OptionWithWarning(LabeledDescriptor):
    def __set__(self, instance, new_option: str):
        if not new_option:  # if new_option is None:
            warnings.warn('Not specifying units is DEPRECATED and will no longer be allowed in the future!',
                          category=DeprecationWarning)
        else:
            instance.__dict__[self.label] = new_option


# ========================================= most important data structures =========================================
class PWscfInput(metaclass=DescriptorOwnerMeta):
    atomic_positions_option = _OptionWithWarning()
    cell_parameters_option = _OptionWithWarning()

    @LazyNamelist
    def control_namelist(self):
        pass

    @LazyNamelist
    def system_namelist(self):
        pass

    @LazyNamelist
    def electrons_namelist(self):
        pass

    @LazyCard
    def cell_parameters(self) -> np.ndarray:
        pass

    @LazyCard
    def atomic_species(self) -> List[NamedTuple]:
        pass

    @LazyCard
    def atomic_positions(self) -> List[NamedTuple]:
        pass

    @LazyWritableProperty
    def k_points(self) -> NamedTuple:
        pass

    @LazyWritableProperty
    def k_points_option(self) -> str:
        pass

    def beautify(self):
        for attr in dir(self):
            if attr.endswith('_namelist'):
                setattr(self, attr, getattr(self, attr).beautify())
        return self

    def to_dict(self):
        d = dict()
        for attr in dir(self):
            if isinstance(getattr(self, attr), NamelistDict):
                d.update({attr: getattr(self, attr).to_dict()})
            elif isinstance(getattr(self, attr), (list, tuple)):
                d.update({attr: getattr(self, attr)})
        return d

    def to_text_file(self, outfile: str = '', path_prefix: Optional[str] = '') -> None:
        outfile_path = path_generator(outfile, path_prefix)

        k_points: KPoints = self.k_points

        with open(outfile_path, 'w') as f:
            f.write("&CONTROL\n")
            f.write(self.control_namelist.to_text())
            f.write("\n/\n&SYSTEM\n")
            f.write(self.system_namelist.to_text())
            f.write("\n/\n&ELECTRONS\n")
            f.write(self.electrons_namelist.to_text())
            f.write("\n/\nCELL_PARAMETERS\n")
            f.write(re.sub("[\[\]]", ' ', np.array2string(self.cell_parameters[0],
                                                          formatter={'float_kind': lambda x: "{:20.10f}".format(x)})))
            f.write("\nATOMIC_SPECIES\n")
            for row in self.atomic_species:
                f.write(' '.join(map(str, row)) + "\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self.atomic_positions_option))
            for row in self.atomic_positions:
                f.write(re.sub("[\[\]]", '', ' '.join(map(str, row)) + "\n"))
            f.write("K_POINTS {{ {0} }}\n".format(self.k_points_option))
            f.write(' '.join(map(str, (k_points.grid + k_points.offsets))))

        print("Object '{0}' is successfully written to file {1}!".format(type(self).__name__,
                                                                         os.path.abspath(outfile_path)))

    def to_json(self) -> None:
        return dumps(self.to_dict())

    def to_json_file(self, outfile: str, path_prefix: Optional[str] = '') -> None:
        outfile_path = path_generator(outfile, path_prefix)
        with open(outfile_path, 'w') as f:
            dump(self.to_dict(), f)
        print("Object '{0}' is successfully written to file {1}!".format(type(self).__name__,
                                                                         os.path.abspath(outfile_path)))


class SCFInput(PWscfInput):
    pass


class VCRelaxInput(PWscfInput):
    @LazyNamelist
    def ions_namelist(self):
        pass

    @LazyNamelist
    def cell_namelist(self):
        pass


class PHononStandardInput:
    @LazyWritableProperty
    def title(self):
        pass

    @LazyNamelist
    def inputph_namelist(self):
        pass

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
            for k, v in self.inputph_namelist.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            if self.single_q_point:
                f.write(' '.join(self.single_q_point))
                f.write("\n")
            elif self.q_points:
                f.write(re.sub("[\[\]]", ' ',
                               np.array2string(self.q_points,
                                               formatter={'float_kind': lambda x: "{:20.10f}".format(x)})))
            else:
                print("Object '{0}' is written to file {1}!".format(self.__name__, os.path.abspath(outfile_path)))

    # TODO: finish this method
    def to_json(self) -> None:
        pass
