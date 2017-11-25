#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from code import InteractiveConsole
from collections import defaultdict

import numpy as np

from basics.pwscf_params import *
from miscellaneous.string import *


class UserFriendlySCFGenerator:
    def __init__(self):
        ic = InteractiveConsole()
        self.raw_parameters = defaultdict(str)
        for k in pseudo_parameters.keys():
            self.raw_parameters[k] = ic.raw_input('Please input {0}: '.format(k))

    def convert_raw_to_standard(self) -> set:
        """
        Convert literally-understandable parameters inputted by the user to Quantum Espresso's standard parameters.

        :return: A set of standard Quantum Espresso parameters.
        """
        official_parameters = set()
        for k, v in self.raw_parameters.items():
            official_parameters.add(self.raw_parameters[k])
        return official_parameters


class SCFGenerator:
    """

    """

    def __init__(self, fixed_in_file: str, variable_in_file: str):
        self.calc_type = 'scf'
        self.fixed_in_file = fixed_in_file
        self.variable_in_file = variable_in_file
        self.parameters = self.read_fixed_input()
        self.control, self.system, self.eletrons, self.ions, self.cell = self._split_parameters()
        self.atomic_species, self.amsmass, self.amspp, self.cell_param, self.atomic_positions = self.read_variable_input()

    def read_fixed_input(self) -> Param:
        """
        This method reads submitters from a given template file. Each line in the file is of form:
            key: value
        Then the method generates a dictionary that contains each different keys and corresponding
        value.

        :return: A dictionary that contains each pair of key and value on each line in the file.
        """
        raw_parameters = defaultdict(str)

        with open(self.fixed_in_file, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                # If a line starts with '#', it will be regarded as a comment,
                # we do not parse this line.
                if stripped_line.startswith('#'):  # If a line starts with '#', continue following code
                    continue
                if not stripped_line:  # If a line is a blank line, jump out of this step and continue next one
                    continue
                # Use ':' as the delimiter, split the line into key and value.
                key, value = stripped_line.split(':', maxsplit=1)
                key: str = key.strip()
                value: str = value.strip()
                if not key:  # If k is ''
                    raise ValueError('The key of {0} is an empty string!'.format(value))
                if not value:  # If v is ''
                    raise ValueError('The value of {0} is an empty string!'.format(key))
                raw_parameters.update({key: value})

        return raw_parameters

    def _split_parameters(self) -> Tuple[DefaultDict[str, List[str]], ...]:
        """
        Split fixed parameters into 5 dictionaries.

        :return: 5 `defaultdict`s.
        """
        # If the keys of `default_parameters` is a subset of that of `self.parameters`
        if default_parameters.keys() <= self.parameters.keys():
            raise KeyError('Your given parameters less than val! The calculation cannot be done!')

        # Repeat n times for a mutable object
        new_control, new_system, new_electrons, new_ions, new_cell = [defaultdict(str) for _ in range(5)]
        # Given a dictionary (`self.parameters`), if a key of it is in `CONTROL`, then append the value
        # of it to `new_control`, etc.
        for k, v in self.parameters.items():
            if k in CONTROL:
                new_control[k] = v
            elif k in SYSTEM:
                new_system[k] = v
            elif k in ELECTRONS:
                new_electrons[k] = v
            elif k in IONS:
                new_ions[k] = v
            elif k in CELL:
                new_cell[k] = v
            else:
                raise KeyError("The key '{0}' you input does not belong to any range!".format(k))
        return new_control, new_system, new_electrons, new_ions, new_cell

    def read_variable_input(self):
        atomic_species = []
        amsmass = []
        amspp = []
        with open(self.variable_in_file, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                # If a line starts with '#', it will be regarded as a comment,
                # we do not parse this line.
                if stripped_line.startswith('#'):  # If a line starts with '#', continue following code
                    continue
                if not stripped_line:  # If a line is a blank line, jump out of this step and continue next one
                    continue
                if 'atomic species' in line:
                    line = f.readline()
                    ams = re.findall("(\w+)\s+(\d+\.?\d+)\s+(\w.*)", line)
                    atomic_species.append(ams[0][0])
                    amsmass.append(ams[0][1])
                    amspp.append(ams[0][2])
                if 'cell parameters' in line:
                    cell_param = np.zeros((3, 3))
                    for i in range(3):
                        sp = f.readline().split()
                        cell_param[i] = strs_to_floats(sp)
                if 'atomic positions' in line:
                    atomic_positions = []
                    line = f.readline()
                    apos = re.findall("(\w+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)", line)
                    atomic_positions.append(apos[0][1:])
        return atomic_species, amsmass, amspp, cell_param, atomic_positions

    def write_scf(self, out_file):
        """
        Generate a PWscf input file for given parameters.

        :param out_file: file that you want to write to as output
        :return:
        """
        with open(out_file, 'w') as f:
            f.write("&CONTROL\n")
            f.write("calculation = '{0}'\n".format(self.calc_type))
            for k, v in self.control.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&SYSTEM\n")
            for k, v in self.system.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&ELECTRONS\n")
            for k, v in self.eletrons.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&IONS\n")
            for k, v in self.ions.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&CELL\n")
            for k, v in self.cell.items():
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")


class VCRelaxGenerator(SCFGenerator):
    def __init__(self, fixed_in_file: str, variable_in_file: str):
        super().__init__(fixed_in_file, variable_in_file)
        self.calc_type = 'vc-relax'
