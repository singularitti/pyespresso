#!/usr/bin/env python3

from typing import *

from f90nml import Parser as NamelistParser

from pyespresso.core.cards import AutomaticKPoints
from pyespresso.core.qe_input import PWInput
from pyespresso.lexers.pwscf import PWscfInputLexer
from pyespresso.meta.text import TextStream

# ========================================= What can be exported? =========================================
__all__ = ['PWInputParser']


class PWInputParser:
    def __init__(self, inp: Optional[str] = None, **kwargs):
        self.__text_stream = TextStream(inp, **kwargs)
        self.lexer = PWscfInputLexer(inp)

    def parse_namelist(self):
        parser = NamelistParser()
        return parser.read(self.__text_stream.stream)

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

    def auto_parse(self) -> PWInput:
        obj = PWInput()
        obj.atomic_species = self.parse_atomic_species()
        obj.atomic_positions = self.parse_atomic_positions()
        obj.k_points = self.parse_k_points()
        obj.cell_parameters = self.parse_cell_parameters()
        obj.eval()
        return obj


class PWOutputParser:
    supported_content = {
        'lattice_parameter': ("lattice parameter \(alat\)\s+=\s*(\d+\.\d+)", float),
        'total_energy': ("!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", float),
        'cell_volume': ("unit-cell volume\s+=\s+(\d+\.\d+)", float),
        'pressure': ("P=\s+(-?\d+\.\d+)", float),
        'kinetic_energy_cutoff': ("kinetic-energy cutoff\s+=\s*(-?\d+\.\d+)", float),
        'charge_density_cutoff': ("charge density cutoff\s+=\s*(-?\d+\.\d+)", float),
        'atoms_num_per_cell': ("number of atoms\/cell\s+=\s*(\d+)", int),
        'atoms_types_num': ("number of atomic INPUTPH_types\s+=\s*(\d+)", int),
        'electrons_num': ("number of electrons\s+=\s*(-?\d+\.\d+)", float),
        'ks_states_num': ("number of Kohn-Sham states\s*=\s*(\d+)", int),
        'mixing_beta': ("mixing beta\s+=\s*(-?\d*\.\d+)", float),
        'nstep': ("nstep\s+=\s*(\d+)", int),
        'iteration_num': ("number of iterations used\s+=\s*(\d+)", int),
        'symmetry_operations_num': ("(\d+)\s+Sym\. Ops.*found", int)
    }
