#!/usr/bin/env python3

from pyque.core.qe_input import PHononStandardInput, PWscfInput
from pyque.lexer.phonon import PHononInputParser
from pyque.lexer.pwscf import PWscfInputParser

# ========================================= What can be exported? =========================================
__all__ = ['build_all', 'PWscfInputBuilder', 'PHononInputBuilder']


def build_all(obj):
    for attr in dir(obj):
        if attr.startswith('build_'):
            getattr(obj, attr)()
    return obj.input_obj


class PWscfInputBuilder:
    def __init__(self, infile):
        self.parser = PWscfInputParser(infile=infile)
        self.input_obj = PWscfInput()

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


class PHononInputBuilder:
    def __init__(self, in_file):
        self.parser = PHononInputParser(in_file)
        self.input_obj = PHononStandardInput()

    def build_inputph_namelist(self):
        self.input_obj.inputph_namelist = self.parser.parse_INPUTPH_namelist()

    def build_title(self):
        self.input_obj.__title__ = self.parser.parse_title()

    def build_single_q_point(self):
        self.input_obj.single_q_point = self.parser.parse_single_q_point()

    def build_q_points(self):
        self.input_obj.q_points = self.parser.parse_q_points()
