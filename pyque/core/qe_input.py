#!/usr/bin/env python3

import os
import re
import warnings
from typing import *
import io

import numpy as np
from json_tricks import dump, dumps
from lazy_property import LazyWritableProperty

from pyque.meta.card import LazyCard
from pyque.meta.descriptors import LabeledDescriptor
from pyque.meta.namelist import LazyNamelist, NamelistDict
from pyque.tools.path_generators import path_generator

# ========================================= What can be exported? =========================================
__all__ = ['is_pwscf_input', 'print_pwscf_input', 'PWscfInput', 'SCFInput', 'VCRelaxInput',
           'PHononInput']

# ========================================= variables declaration =========================================
_input_attr = {'CONTROL': 'control_namelist', 'SYSTEM': 'system_namelist', 'ELECTRONS': 'electrons_namelist'}


# ========================================= define useful functions =========================================
def is_pwscf_input(obj: object):
    if isinstance(obj, PWscfInput):
        return True
    return False


def print_pwscf_input(obj: object):
    if is_pwscf_input(obj):
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
class PWscfInput:
    @LazyNamelist
    def control_namelist(self) -> NamelistDict:
        """
        Input variables that control the flux of the calculation and the amount of I/O on disk and on the screen.

        :return:
        """
        return NamelistDict()

    @LazyNamelist
    def system_namelist(self) -> NamelistDict:
        """
        Input variables that specify the system under study.

        :return:
        """
        return NamelistDict()

    @LazyNamelist
    def electrons_namelist(self) -> NamelistDict:
        """
        input variables that control the algorithms used to reach the self-consistent solution of KS
        equations for the electrons.

        :return:
        """
        return NamelistDict()

    @LazyNamelist
    def ions_namelist(self) -> NamelistDict:
        return NamelistDict()

    @LazyNamelist
    def cell_namelist(self) -> NamelistDict:
        return NamelistDict()

    @LazyCard
    def cell_parameters(self) -> Dict:
        return {'option': None, 'value': None}

    @LazyCard
    def atomic_species(self) -> Dict:
        """
        name, mass and pseudopotential used for each atomic species present in the system
        :return:
        """
        return {'option': None, 'value': None}

    @LazyCard
    def atomic_positions(self) -> Dict:
        """
        type and coordinates of each atom in the unit cell
        :return:
        """
        return {'option': None, 'value': None}

    @LazyWritableProperty
    def k_points(self) -> Dict:
        """
        coordinates and weights of the k-points used for BZ integration
        :return:
        """
        return {'option': None, 'value': None}

    def eval(self):
        for attr in dir(self):
            if hasattr(attr, 'eval'):
                if getattr(self, attr) is None:
                    continue
                getattr(self, '_' + attr).eval()

    def to_dict(self):
        d = dict()
        for attr in dir(self):
            if isinstance(getattr(self, attr), NamelistDict):
                d.update({attr: getattr(self, attr).to_dict()})
            elif isinstance(getattr(self, attr), (list, tuple)):
                d.update({attr: getattr(self, attr)})
        return d

    def to_text(self):
        s = list()
        s.append("&CONTROL")
        s.append(self.control_namelist.to_text())
        s.append("/\n&SYSTEM")
        s.append(self.system_namelist.to_text())
        s.append("/\n&ELECTRONS")
        s.append(self.electrons_namelist.to_text())
        s.append("/\n&IONS")
        if self.ions_namelist is None:
            pass
        else:
            s.append(self.ions_namelist.to_text())
        s.append("/\n&CELL")
        if self.cell_namelist is None:
            pass
        else:
            s.append(self.cell_namelist.to_text())
        s.append("/\nATOMIC_SPECIES")
        for row in self.atomic_species['value']:
            s.append(str(row))
        s.append("ATOMIC_POSITIONS {{ {0} }}".format(self.atomic_positions['option']))
        for row in self.atomic_positions['value']:
            s.append(str(row))
        s.append("K_POINTS {{ {0} }}".format(self.k_points['option']))
        s.append(str(self.k_points['value']))
        if self.cell_parameters is None:
            pass
        else:
            s.append("CELL_PARAMETERS {{ {0} }}".format(self.cell_parameters['option']))
            s.append(re.sub("[\[\]']", ' ', np.array2string(self.cell_parameters['value'],
                                                            formatter={'float_kind': lambda x: '{:20.10f}'.format(x)})))
        return '\n'.join(s)

    def to_text_file(self, outfile: str = '', path_prefix: Optional[str] = '') -> None:
        outfile_path = path_generator(outfile, path_prefix)
        with open(outfile_path, 'w') as f:
            f.write(self.to_text())
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
        """
        needed when ATOMS MOVE! IGNORED otherwise ! input variables that control ionic motion in
        molecular dynamics run or structural relaxation.

        :return:
        """
        pass

    @LazyNamelist
    def cell_namelist(self):
        """
        needed when CELL MOVES! IGNORED otherwise ! input variables that control the cell-shape evolution
        in a variable-cell-shape MD or structural relaxation.

        :return:
        """
        pass


class PHononInput:
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
