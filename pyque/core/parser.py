#!/usr/bin/env python3

from pyque.core.qe_input import PHononInput, PWscfInput
from pyque.lexer.pwscf import PWscfInputLexer

# ========================================= What can be exported? =========================================
__all__ = ['auto_build', 'PWscfInputParser', 'PHononInputParser']


def auto_build(obj):
    if isinstance(obj, PWscfInputParser):
        obj.auto_build()
    return obj.input_obj.eval()


class PWscfInputParser:
    def __init__(self, infile):
        self.lexer = PWscfInputLexer(infile=infile)
        self.input_obj = PWscfInput()
        self.parse_dict = {
            '&CONTROL': 'parse_control_namelist',
            '&SYSTEM': 'parse_system_namelist',
            '&ELECTRONS': 'parse_electrons_namelist',
            'ATOMIC_SPECIES': 'parse_atomic_species',
            'ATOMIC_POSITIONS': 'parse_atomic_positions',
            'K_POINTS': 'parse_k_points',
            'CELL_PARAMETERS': 'parse_cell_parameters'
        }

    def build_control_namelist(self):
        self.input_obj.control_namelist = self.lexer.parse_control_namelist()

    def build_system_namelist(self):
        self.input_obj.system_namelist = self.lexer.parse_system_namelist()

    def build_electrons_namelist(self):
        self.input_obj.electrons_namelist = self.lexer.parse_electrons_namelist()

    def build_atomic_species(self):
        self.input_obj.atomic_species = self.lexer.parse_atomic_species()

    def build_atomic_positions(self):
        self.input_obj.atomic_positions, self.input_obj.atomic_positions_option = self.lexer.parse_atomic_positions()

    def build_k_points(self):
        self.input_obj.k_points, self.input_obj.k_points_option = self.lexer.parse_k_points()

    def build_cell_parameters(self):
        self.input_obj.cell_parameters = self.lexer.parse_cell_parameters()

    def auto_build(self):
        provided_namelists = self.lexer.namelist_identifier_positions.keys()
        provided_cards = self.lexer.card_identifier_positions.keys()
        methods_involved = [self.parse_dict[identifiers] for identifiers in provided_namelists] + [
            self.parse_dict[identifiers] for identifiers in provided_cards]
        for method in methods_involved:
            pass


class PHononInputParser:
    def __init__(self, in_file):
        self.parser = PHononInputParser(in_file)
        self.input_obj = PHononInput()

    def build_inputph_namelist(self):
        self.input_obj.inputph_namelist = self.parser.parse_inputph_namelist()

    def build_title(self):
        self.input_obj.__title__ = self.parser.parse_title()

    def build_single_q_point(self):
        self.input_obj.single_q_point = self.parser.parse_single_q_point()

    def build_q_points(self):
        self.input_obj.q_points = self.parser.parse_q_points()
