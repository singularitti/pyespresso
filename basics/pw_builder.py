#!/usr/bin/env python3
# created at Dec 7, 2017 1:46 AM by Qi Zhang

from basics.pw_params import *
from basics.pwscf import PWStandardInput
from miscellaneous.string import *
from readers.pwscf import PWInputParser


def build_pw_input(in_file: str) -> PWStandardInput:
    pwib = PWInputBuilder(in_file)
    pwib.build_all()
    return pwib.input_obj


class PWInputBuilder:
    def __init__(self, in_file):
        self.parser = PWInputParser(in_file)
        self.input_obj = PWStandardInput()

    def build_CONTROL(self):
        self.input_obj.CONTROL = self.parser.parse_CONTROL_namelist()

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
        CONTROL_default_parameters = CONTROL_namelist.default_parameters
        for k, v in self.pw.CONTROL.items():
            if CONTROL_default_parameters[k][1] is float:
                d.update({k: str_to_float(v)})
            else:
                d.update({k: CONTROL_default_parameters[k][1](v)})
        return d

    def fancy_SYSTEM(self):
        d = {}
        SYSTEM_default_parameters = SYSTEM_namelist.default_parameters
        for k, v in self.pw.SYSTEM.items():
            if '(' in k:
                # Only take the part before the first '('
                k_prefix = re.match("(\w+)\(?(\d*)\)?", k, flags=re.IGNORECASE).group(1)
            else:
                k_prefix = k
            print(k_prefix)
            if SYSTEM_default_parameters[k_prefix][1] is float:
                d.update({k: str_to_float(v)})
            else:
                d.update({k: SYSTEM_default_parameters[k_prefix][1](v)})
        return d

    def fancy_ELECTRONS(self):
        d = {}
        ELECTRONS_default_parameters = ELECTRONS_namelist.default_parameters
        for k, v in self.pw.ELECTRONS.items():
            if ELECTRONS_default_parameters[k][1] is float:
                d.update({k: str_to_float(v)})
            else:
                d.update({k: ELECTRONSParameter(k, v)})
        return d

    def fancy_input(self):
        self.pw.CONTROL = self.fancy_CONTROL()
        self.pw.SYSTEM = self.fancy_SYSTEM()
        self.pw.ELECTRONS = self.fancy_ELECTRONS()
        return self.pw
