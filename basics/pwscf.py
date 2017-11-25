#!/usr/bin/env python3
# created at Nov 21, 2017 2:13 AM by Qi Zhang

import re
from collections import namedtuple
from typing import *

import numpy as np

from basics.tree import Tree
from miscellaneous.dictionary import merge_dicts

k_mesh = namedtuple('k_mesh', ['grid', 'shift', 'option'])


class SCFStandardInput(Tree):
    def __init__(self, in_file):
        super().__init__()
        self.in_file = in_file
        self.control_card = {}
        self.system_card = {}
        self.electrons_card = {}
        self.cell_parameters: Dict[str, np.ndarray] = {'CELL_PARAMETERS': np.zeros((3, 3))}
        self.k_mesh: Dict[str, NamedTuple] = {'K_POINTS': k_mesh(grid=[1, 1, 1], shift=[0, 0, 0], option='')}
        self.atomic_species = {'ATOMIC_SPECIES': []}
        self.atomic_positions = {'ATOMIC_POSITIONS': [], 'option': ''}

    @property
    def flattened_tree(self) -> Dict[str, dict]:
        return merge_dicts(
            self.control_card, self.system_card, self.electrons_card, self.cell_parameters, self.k_mesh
        )

    def write_to_file(self, out_file: str):
        k_mesh = self.k_mesh['K_POINTS']
        with open(out_file, 'w') as f:
            f.write("&CONTROL\n")
            for k, v in self.control_card.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&SYSTEM\n")
            for k, v in self.system_card.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n&ELECTRONS\n")
            for k, v in self.electrons_card.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\nCELL_PARAMETERS\n")
            f.write(
                re.sub("[\[\]]", ' ', np.array2string(self.cell_parameters['CELL_PARAMETERS'],
                                                      formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
            )
            f.write("\nATOMIC_SPECIES\n")
            for row in self.atomic_species['ATOMIC_SPECIES']:
                for i in row:
                    f.write("{0}   ".format(i))
                f.write("\n")
            f.write("ATOMIC_POSITIONS {{ {0} }}\n".format(self.atomic_positions['option']))
            for row in self.atomic_positions['ATOMIC_POSITIONS']:
                for i in row:
                    f.write("{0}   ".format(i))
                f.write("\n")
            f.write("K_POINTS {{ {0} }}\n".format(k_mesh.option))
            f.write(' '.join(map(str, (k_mesh.grid + k_mesh.shift))))
            f.write("\n")

    def __repr__(self):
        from beeprint import pp
        return pp(self.flattened_tree, output=False)

# class PWscfInteractiveInput(InteractiveConsole):
#     def aaa(self):
#         psi = PWscfStandardInput()
#         k = self.raw_input('Please input an option name: ')
#         v = self.raw_input('Please give its value: ')
#         if k in default_parameters:
#
