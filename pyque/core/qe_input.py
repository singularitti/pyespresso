#!/usr/bin/env python3

import os
import re
import warnings
from collections import namedtuple
from typing import *

import numpy as np
from lazy_property import LazyWritableProperty

from pyque.meta.descriptors import LabeledDescriptor, DescriptorOwnerMeta
from pyque.meta.namelist import ELECTRONS_NAMELIST, CONTROL_NAMELIST, SYSTEM_NAMELIST, DefaultParameters
from pyque.meta.parameter import to_qe_str
from pyque.util.path_generators import path_generator
from pyque.util.strings import *
from pyque.lexer.simple import ValueWithComment

# ========================================= What can be exported? =========================================
__all__ = ['is_pw_input', 'print_pw_input', 'PWscfStandardInput', 'SCFStandardInput', 'VCRelaxStandardInput',
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
control_default_parameters: DefaultParameters = CONTROL_NAMELIST.default_parameters
system_default_parameters: DefaultParameters = SYSTEM_NAMELIST.default_parameters
electrons_default_parameters: DefaultParameters = ELECTRONS_NAMELIST.default_parameters


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


class _CONTROLNamelist(LabeledDescriptor):
    pass


# ========================================= most important data structures =========================================
class PWscfStandardInput(metaclass=DescriptorOwnerMeta):
    atomic_positions_option = _OptionWithWarning()
    cell_parameters_option = _OptionWithWarning()

    @LazyWritableProperty
    def control_namelist(self):
        pass

    @LazyWritableProperty
    def system_namelist(self):
        pass

    @LazyWritableProperty
    def electrons_namelist(self):
        pass

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
            for k, v in self.control_namelist.items():
                f.write("{0} = {1}\n".format(k, to_qe_str(v)))
            f.write("/\n&SYSTEM\n")
            for k, v in self.system_namelist.items():
                f.write("{0} = {1}\n".format(k, to_qe_str(v)))
            f.write("/\n&ELECTRONS\n")
            for k, v in self.electrons_namelist.items():
                f.write("{0} = {1}\n".format(k, to_qe_str(v)))
            f.write("/\nCELL_PARAMETERS\n")
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

        print("Object '{0}' is written to file {1}!".format(type(self).__name__, os.path.abspath(outfile_path)))

    # TODO: finish this method
    def to_json(self) -> None:
        pass

    def beautify(self):
        d = {}
        for k, v in self.control_namelist.items():
            if control_default_parameters[k][1] is float:
                d.update({k: ValueWithComment(string_to_general_float(v.value), v.comment)})
            else:
                d.update({k: ValueWithComment(control_default_parameters[k][1](v.value), v.comment)})
        self.control_namelist = d
        d = {}
        for k, v in self.system_namelist.items():
            if '(' in k:
                # Only take the part before the first '(' to deal with names like 'celldm(1)'.
                k_prefix = re.match("(\w+)\(?(\d*)\)?", k, flags=re.IGNORECASE).group(1)
            else:
                k_prefix = k
            if system_default_parameters[k_prefix][1] is float:
                d.update({k: ValueWithComment(string_to_general_float(v.value), v.comment)})
            else:
                d.update({k: ValueWithComment(system_default_parameters[k_prefix][1](v.value), v.comment)})
        self.system_namelist = d
        d = {}
        for k, v in self.electrons_namelist.items():
            if electrons_default_parameters[k][1] is float:
                d.update({k: ValueWithComment(string_to_general_float(v.value), v.comment)})
            else:
                d.update({k: ValueWithComment(electrons_default_parameters[k][1](v.value), v.comment)})
        self.electrons_namelist = d
        return self


class SCFStandardInput(PWscfStandardInput):
    pass


class VCRelaxStandardInput(PWscfStandardInput):
    @LazyWritableProperty
    def ions_namelist(self):
        pass

    @LazyWritableProperty
    def cell_namelist(self):
        pass

    # cell_namelist.__name__ = 'CELL Namelist'


class PHononStandardInput:
    @LazyWritableProperty
    def title(self):
        pass

    @LazyWritableProperty
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
