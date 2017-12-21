#!/usr/bin/env python3
# created at Dec 7, 2017 1:46 AM by Qi Zhang

from basics.pwscf import PWStandardInput
from readers.pwscf import PWInputParser
from basics.pw_params import *
from miscellaneous.string import *


def build_pw_input(in_file: str) -> PWStandardInput:
    pwib = PWInputBuilder(in_file)
    pwib.build_all()
    return pwib.input_obj


class PWInputBuilder:
    def __init__(self, in_file):
        self.parser = PWInputParser(in_file)
        self.input_obj = PWStandardInput()

    def build_CONTROL(self):
        self.input_obj.CONTROL = self.parser.parse_control_namelist()

    def build_SYSTEM(self):
        self.input_obj.SYSTEM = self.parser.parse_system_namelist()

    def build_ELECTRONS(self):
        self.input_obj.ELECTRONS = self.parser.parse_electrons_namelist()

    def build_CELL_PARAMETERS(self):
        self.input_obj.CELL_PARAMETERS = self.parser.parse_cell_parameters()

    def build_K_POINTS(self):
        self.input_obj.K_POINTS, self.input_obj.K_POINTS_option = self.parser.parse_k_points()

    def build_ATOMIC_SPECIES(self):
        self.input_obj.ATOMIC_SPECIES = self.parser.parse_atomic_species()

    def build_ATOMIC_POSITIONS(self):
        self.input_obj.ATOMIC_POSITIONS, self.input_obj.ATOMIC_POSITIONS_option = self.parser.parse_atomic_positions()

    def build_all(self):
        self.build_CONTROL()
        self.build_SYSTEM()
        self.build_ELECTRONS()
        self.build_CELL_PARAMETERS()
        self.build_K_POINTS()
        self.build_ATOMIC_SPECIES()
        self.build_ATOMIC_POSITIONS()


class PWInputFancier:
    def __init__(self, obj: PWStandardInput):
        self.pw = obj

    def fancy_CONTROL(self):
        d = {}
        for k, v in self.pw.CONTROL.items():
            if CONTROL_values[k][1] is float:
                d.update({k: str_to_float(v)})
            else:
                d.update({k: CONTROL_values[k][1](v)})
        return d

    def fancy_ELECTRONS(self):
        d = {}
        for k, v in self.pw.ELECTRONS.items():
            if CONTROL_values[k][1] is float:
                d.update({k: str_to_float(v)})
            else:
                d.update({k: CONTROL_values[k][1](v)})
        return d
