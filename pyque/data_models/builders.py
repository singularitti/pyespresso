#!/usr/bin/env python3

from pyque.data_models.qe_input import PHononStandardInput, PWscfStandardInput
from pyque.meta.namelist import ELECTRONS_NAMELIST, CONTROL_NAMELIST, SYSTEM_NAMELIST
from pyque.meta.parameter import ELECTRONSParameter
from pyque.miscellaneous.strings import *
from pyque.parsers.phonon import PHononInputParser
from pyque.parsers.pwscf import PWscfStandardInputParser

# ========================================= What can be exported? =========================================
__all__ = ['build_pw_input', 'PWscfStandardInputBuilder', 'PWInputFancier', 'PHononInputBuilder']


def build_pw_input(in_file: str) -> PWscfStandardInput:
    pwib = PWscfStandardInputBuilder(in_file)
    pwib.build_all()
    return pwib.input_obj


class PWscfStandardInputBuilder:
    def __init__(self, infile):
        self.parser = PWscfStandardInputParser(infile=infile)
        self.input_obj = PWscfStandardInput()

    def build_CONTROL(self):
        self.input_obj.control_namelist = self.parser.parse_control_namelist()

    def build_SYSTEM(self):
        self.input_obj.system_namelist = self.parser.parse_system_namelist()

    def build_ELECTRONS(self):
        self.input_obj.electrons_namelist = self.parser.parse_electrons_namelist()

    def build_CELL_PARAMETERS(self):
        self.input_obj.cell_parameters = self.parser.parse_cell_parameters()

    def build_K_POINTS(self):
        self.input_obj.k_points, self.input_obj.k_points_option = self.parser.parse_k_points()

    def build_ATOMIC_SPECIES(self):
        self.input_obj.atomic_species = self.parser.parse_atomic_species()

    def build_ATOMIC_POSITIONS(self):
        self.input_obj.atomic_positions, self.input_obj.atomic_positions_option = self.parser.parse_atomic_positions()

    def build_all(self):
        self.build_CONTROL()
        self.build_SYSTEM()
        self.build_ELECTRONS()
        self.build_CELL_PARAMETERS()
        self.build_K_POINTS()
        self.build_ATOMIC_SPECIES()
        self.build_ATOMIC_POSITIONS()


class PWInputFancier:
    def __init__(self, obj: PWscfStandardInput):
        self.pw = obj

    def fancy_CONTROL(self):
        d = {}
        CONTROL_default_parameters = CONTROL_NAMELIST.default_parameters
        for k, v in self.pw.control_namelist.items():
            if CONTROL_default_parameters[k][1] is float:
                d.update({k: str_to_float(v.value)})
            else:
                d.update({k: CONTROL_default_parameters[k][1](v.value)})
        return d

    def fancy_SYSTEM(self):
        d = {}
        SYSTEM_default_parameters = SYSTEM_NAMELIST.default_parameters
        for k, v in self.pw.system_namelist.items():
            if '(' in k:
                # Only take the part before the first '('
                k_prefix = re.match("(\w+)\(?(\d*)\)?", k, flags=re.IGNORECASE).group(1)
            else:
                k_prefix = k
            print(k_prefix)
            if SYSTEM_default_parameters[k_prefix][1] is float:
                d.update({k: str_to_float(v.value)})
            else:
                d.update({k: SYSTEM_default_parameters[k_prefix][1](v.value)})
        return d

    def fancy_ELECTRONS(self):
        d = {}
        ELECTRONS_default_parameters = ELECTRONS_NAMELIST.default_parameters
        for k, v in self.pw.electrons_namelist.items():
            if ELECTRONS_default_parameters[k][1] is float:
                d.update({k: str_to_float(v.value)})
            else:
                d.update({k: ELECTRONSParameter(k, v.value)})
        return d

    def fancy_input(self):
        self.pw.control_namelist = self.fancy_CONTROL()
        self.pw.system_namelist = self.fancy_SYSTEM()
        self.pw.electrons_namelist = self.fancy_ELECTRONS()
        return self.pw


class PHononInputBuilder:
    def __init__(self, in_file):
        self.parser = PHononInputParser(in_file)
        self.input_obj = PHononStandardInput()

    def build_INPUTPH(self):
        self.input_obj.inputph_namelist = self.parser.parse_INPUTPH_namelist()

    def build_title(self):
        self.input_obj.__title__ = self.parser.parse_title()

    def build_single_q_point(self):
        self.input_obj.single_q_point = self.parser.parse_single_q_point()

    def build_q_points(self):
        self.input_obj.q_points = self.parser.parse_q_points()

    def build_all(self):
        self.build_INPUTPH()
        # self.build_title()
        self.build_single_q_point()
        self.build_q_points()
