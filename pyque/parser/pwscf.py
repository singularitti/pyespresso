#!/usr/bin/env python3

from typing import *

from pyque.core.cards import AutomaticKPoints
from pyque.core.qe_input import PWscfInput
from pyque.lexer.pwscf import PWscfInputLexer
from pyque.meta.namelist import *

# ========================================= What can be exported? =========================================
__all__ = ['PWscfInputParser']


class PWscfInputParser:
    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None):
        self.lexer = PWscfInputLexer(instream, infile)

    def parse_namelist(self, group_name: str):
        return {'CONTROL': self.parse_control_namelist, 'SYSTEM': self.parse_system_namelist,
                'ELECTRONS': self.parse_electrons_namelist,
                'IONS': self.parse_ions_namelist, 'CELL': self.parse_cell_namelist}[group_name]()

    def parse_control_namelist(self):
        d = self.lexer.lex_namelist('CONTROL')
        return NamelistDict({k: CONTROLNamelistVariable(k, v) for k, v in d.items()})

    def parse_system_namelist(self):
        d = self.lexer.lex_namelist('SYSTEM')
        return NamelistDict({k: SYSTEMNamelistVariable(k, v) for k, v in d.items()})

    def parse_electrons_namelist(self):
        d = self.lexer.lex_namelist('ELECTRONS')
        return NamelistDict({k: ELECTRONSNamelistVariable(k, v) for k, v in d.items()})

    def parse_ions_namelist(self):
        d = self.lexer.lex_namelist('IONS')
        return NamelistDict({k: IONSNamelistVariable(k, v) for k, v in d.items()})

    def parse_cell_namelist(self):
        d = self.lexer.lex_namelist('CELL')
        return NamelistDict({k: CELLNamelistVariable(k, v) for k, v in d.items()})

    def parse_card(self, card_name: str):
        return {'ATOMIC_SPECIES': self.parse_atomic_species, 'ATOMIC_POSITIONS': self.parse_atomic_positions,
                'K_POINTS': self.parse_k_points, 'CELL_PARAMETERS': self.parse_cell_parameters,
                'OCCUPATIONS': NotImplemented, 'CONSTRAINTS': NotImplemented,
                'ATOMIC_FORCES': NotImplemented}[card_name]

    def parse_atomic_species(self):
        return [_.eval() for _ in self.lexer.lex_atomic_species()]

    def parse_atomic_positions(self):
        atomic_positions, atomic_positions_option = self.lexer.lex_atomic_positions()
        return [_.eval() for _ in atomic_positions], atomic_positions_option

    def parse_k_points(self):
        k_points, k_points_option = self.lexer.lex_k_points()
        if isinstance(k_points, AutomaticKPoints):
            k_points = k_points.eval()
        return k_points, k_points_option

    def parse_cell_parameters(self):
        return self.lexer.lex_cell_parameters()

    def auto_parse(self):
        namelists_found: Optional[Set[str]] = self.lexer.namelists_found
        cards_found: Optional[Set[str]] = self.lexer.cards_found
        obj = PWscfInput()
        for _ in namelists_found:
            obj._receive(self.parse_namelist(_))
        obj.atomic_species = self.parse_atomic_species()
        obj.atomic_positions = self.parse_atomic_positions()
        obj.k_points = self.parse_k_points()
        for _ in (cards_found - {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'}):
            obj._receive(self.parse_card(_))
