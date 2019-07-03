#!/usr/bin/env python3

import os
import re
from typing import *

import addict
import numpy as np
from json_tricks import dump, dumps
from lazy_property import LazyWritableProperty

from pyespresso.core.namelist import LazyNamelist, NamelistDict
from pyespresso.tools.os_path import path_generator

# ========================================= What can be exported? =========================================
__all__ = ['is_pw_input', 'print_pw_input', 'PWInput', 'PHononInput']

# ========================================= variables declaration =========================================
_input_attr = {'CONTROL': 'control_namelist', 'SYSTEM': 'system_namelist', 'ELECTRONS': 'electrons_namelist'}


# ========================================= define useful functions =========================================
def is_pw_input(obj: object):
    if isinstance(obj, PWInput):
        return True
    return False


def print_pw_input(obj: object):
    if is_pw_input(obj):
        try:
            from beeprint import pp
            pp(obj)
        except ModuleNotFoundError:
            print(obj)
    else:
        raise TypeError("{0} is not a {1} object!".format(obj, obj.__class__.__name__))


# ========================================= define useful data structures =========================================


# ========================================= most important data structures =========================================
class PWInput(addict.Dict):
    __slots__ = ('namelists', 'cards', 'eval', 'to_dict', 'to_text', 'to_text_file', 'to_json', 'to_json_file')

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
