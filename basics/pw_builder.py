#!/usr/bin/env python3
# created at Dec 7, 2017 1:46 AM by Qi Zhang

from basics.pwscf import PWStandardInput
from readers.pwscf import PWInputParser


def build_pw_input(in_file: str) -> PWStandardInput:
    pwib = PWInputBuilder(in_file)
    pwib.build_all()
    return pwib.input_obj


class PWInputBuilder:
    def __init__(self, in_file):
        self.parser = PWInputParser(in_file)
        self.input_obj = PWStandardInput()

    def build_control_namelist(self):
        self.input_obj.control_namelist = self.parser.parse_control_namelist()

    def build_system_namelist(self):
        self.input_obj.system_namelist = self.parser.parse_system_namelist()

    def build_electrons_namelist(self):
        self.input_obj.electrons_namelist = self.parser.parse_electrons_namelist()

    def build_cell_parameters(self):
        self.input_obj.cell_parameters = self.parser.parse_cell_parameters()

    def build_k_points(self):
        self.input_obj.k_points, self.input_obj.k_points_option = self.parser.parse_k_points()

    def build_atomic_species(self):
        self.input_obj.atomic_species = self.parser.parse_atomic_species()

    def build_atomic_positions(self):
        self.input_obj.atomic_positions, self.input_obj.atomic_positions_option = self.parser.parse_atomic_positions()

    def build_all(self):
        self.build_control_namelist()
        self.build_system_namelist()
        self.build_electrons_namelist()
        self.build_cell_parameters()
        self.build_k_points()
        self.build_atomic_species()
        self.build_atomic_positions()
