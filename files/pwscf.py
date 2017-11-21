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
        super().__init__({})
        self.in_file = in_file
        self._control_card = {}
        self._system_card = {}
        self._electrons_card = {}
        self._cell_parameters: np.ndarray = np.zeros((3, 3))
        self._k_mesh: NamedTuple = k_mesh(grid=(1, 1, 1), shift=(0, 0, 0))

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
    def k_mesh(self) -> NamedTuple:
        return self._k_mesh

    @k_mesh.setter
    def k_mesh(self, val: NamedTuple):
        self.k_mesh = val

    @property
    def cell_parameters(self) -> np.ndarray:
        return self._cell_parameters

    @cell_parameters.setter
    def cell_parameters(self, val: np.ndarray):
        self._cell_parameters = val

    @property
    def flattened_tree(self) -> Dict:
        print(self._cell_parameters)
        return merge_dicts(self._control_card, self._system_card, self._electrons_card
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
