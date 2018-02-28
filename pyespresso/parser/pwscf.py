#!/usr/bin/env python3

from typing import *

from pyespresso.core.cards import AutomaticKPoints
from pyespresso.core.qe_input import PWscfInput
from pyespresso.lexer.pwscf import PWscfInputLexer
from pyespresso.core.namelist import namelist_variable, NamelistDict

# ========================================= What can be exported? =========================================
__all__ = ['PWscfInputParser']


class PWscfInputParser:
    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None):
        self.lexer = PWscfInputLexer(instream, infile)

    def parse_namelist(self, group_name: str):
        d = self.lexer.lex_namelist(group_name)
        # You do not need to worry about ``NamelistVariable`` here if you use ``==``.
        return NamelistDict({k: namelist_variable(group_name, k, v) for k, v in d.items()})

    def parse_control_namelist(self):
        return self.parse_namelist('CONTROL')

    def parse_system_namelist(self):
        return self.parse_namelist('SYSTEM')

    def parse_electrons_namelist(self):
        return self.parse_namelist('ELECTRONS')

    def parse_ions_namelist(self):
        return self.parse_namelist('IONS')

    def parse_cell_namelist(self):
        return self.parse_namelist('CELL')

    def parse_card(self, card_name: str):
        return {'ATOMIC_SPECIES': self.parse_atomic_species, 'ATOMIC_POSITIONS': self.parse_atomic_positions,
                'K_POINTS': self.parse_k_points, 'CELL_PARAMETERS': self.parse_cell_parameters,
                'OCCUPATIONS': NotImplemented, 'CONSTRAINTS': NotImplemented,
                'ATOMIC_FORCES': NotImplemented}[card_name]

    def parse_atomic_species(self):
        return {'value': [_.eval() for _ in self.lexer.lex_atomic_species()]}

    def parse_atomic_positions(self):
        atomic_positions, atomic_positions_option = self.lexer.lex_atomic_positions()
        return {'value': [_.eval() for _ in atomic_positions], 'option': atomic_positions_option}

    def parse_k_points(self):
        k_points, k_points_option = self.lexer.lex_k_points()
        if isinstance(k_points, AutomaticKPoints):
            k_points = k_points.eval()
        return {'value': k_points, 'option': k_points_option}

    def parse_cell_parameters(self):
        if self.lexer.lex_cell_parameters() is None:
            return None
        cell_params, option = self.lexer.lex_cell_parameters()
        return {'value': cell_params, 'option': option}

    def auto_parse(self) -> PWscfInput:
        namelists_found: Optional[Set[str]] = self.lexer.namelists_found
        cards_found: Optional[Set[str]] = self.lexer.cards_found
        obj = PWscfInput()
        for _ in namelists_found:
            setattr(obj, _.lower() + '_namelist', self.parse_namelist(_))
        obj.atomic_species = self.parse_atomic_species()
        obj.atomic_positions = self.parse_atomic_positions()
        obj.k_points = self.parse_k_points()
        obj.cell_parameters = self.parse_cell_parameters()
        # for _ in (cards_found - {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'}):
        #     setattr(obj, _.lower(), self.parse_card(_))
        obj.eval()
        return obj
