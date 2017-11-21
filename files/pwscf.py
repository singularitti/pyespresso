#!/usr/bin/env python3
# created at Nov 21, 2017 2:13 AM by Qi Zhang

from collections import namedtuple
from typing import *

import numpy as np

from miscellaneous.dictionary import merge_dicts
from miscellaneous.tree import Tree

k_mesh = namedtuple('k_mesh', ['grid', 'shift'])


class PWscfStandardInput(Tree):
    def __init__(self, in_file):
        super().__init__()
        self.in_file = in_file
        self._control_card = {}
        self._system_card = {}
        self._electrons_card = {}
        self._cell_parameters: Dict[str, np.ndarray] = {'CELL_PARAMETERS': np.zeros((3, 3))}
        self._k_mesh: Dict[str, NamedTuple] = {'K_POINTS': k_mesh(grid=(1, 1, 1), shift=(0, 0, 0))}

    @property
    def control_card(self) -> Dict:
        return self._control_card

    @control_card.setter
    def control_card(self, val: Dict):
        self._control_card = val

    @property
    def system_card(self) -> Dict:
        return self._system_card

    @system_card.setter
    def system_card(self, val: Dict):
        self._system_card = val

    @property
    def electrons_card(self) -> Dict:
        return self._electrons_card

    @electrons_card.setter
    def electrons_card(self, val: Dict):
        self._electrons_card = val

    @property
    def k_mesh(self) -> Dict[str, NamedTuple]:
        return self._k_mesh

    @k_mesh.setter
    def k_mesh(self, val: NamedTuple):
        self._k_mesh = val

    @property
    def cell_parameters(self) -> Dict[str, np.ndarray]:
        return self._cell_parameters

    @cell_parameters.setter
    def cell_parameters(self, val: np.ndarray):
        self._cell_parameters = val

    @property
    def flattened_tree(self) -> Dict[str, dict]:
        return merge_dicts(
            self._control_card, self._system_card, self._electrons_card, self._cell_parameters, self._k_mesh
        )

    def __str__(self):
        return str(self.flattened_tree)

    __repr__ = __str__

# class PWscfInteractiveInput(InteractiveConsole):
#     def aaa(self):
#         psi = PWscfStandardInput()
#         k = self.raw_input('Please input an option name: ')
#         v = self.raw_input('Please give its value: ')
#         if k in default_parameters:
#
