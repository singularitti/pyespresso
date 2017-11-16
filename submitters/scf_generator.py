#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from code import InteractiveConsole
from collections import defaultdict
from typing import *

from submitters.pw_paramlist import *

# Type aliasing
Param = DefaultDict[str, List[str]]


class SCFGenerator:
    def __init__(self, in_file: str, calc_type: str):
        self.in_file = in_file
        self.calc_type = calc_type
        self.parameters = self.read_pwscf_input()
        self.control, self.system, self.eletrons, self.ions, self.cell = self._split_parameters()

    def read_pwscf_input(self) -> Param:
        """
        This method reads submitters from a given template file. Each line in the file is of form:
            key: value
        Then the method generates a dictionary that contains each different keys and corresponding
        value.

        :return: A dictionary that contains each pair of key and value on each line in the file.
        """
        parameters = defaultdict(list)

        with open(self.in_file, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                # If a line starts with '#', it will be regarded as a comment,
                # we do not parse this line.
                if not stripped_line.startswith('#'):
                    continue
                if stripped_line:  # If this line is not a blank line
                    continue
                # Use ':' as the delimiter, split the line into key and value.
                key, value = stripped_line.split(':', maxsplit=1)
                key: str = key.strip()
                value: str = value.strip()
                if not key:  # If k is ''
                    raise ValueError('The key of {0} is an empty string!'.format(value))
                if not value:  # If v is ''
                    raise ValueError('The value of {0} is an empty string!'.format(key))
                parameters.update({key: value})
                print(parameters)

        return parameters

    def _split_parameters(self) -> Tuple[Param, Param, Param, Param, Param]:
        """


        :return:
        """
        # If the keys of `default_parameters` is a subset of that of `self.parameters`
        if default_parameters.keys() <= self.parameters.keys():
            raise KeyError('Your given parameters less than default! The calculation cannot be done!')
        new_control, new_system, new_electrons, new_ions, new_cell = [defaultdict(list)] * 5
        # Given a dictionary (`self.parameters`), if a key of it is in `CONTROL`, then append the value
        # of it to `new_control`, etc.
        for k, v in self.parameters:
            if k in CONTROL:
                print(1)
                new_control[k].append(v)
            if k in SYSTEM:
                new_system[k].append(v)
            if k in ELECTRONS:
                new_electrons[k].append(v)
            if k in IONS:
                new_ions[k].append(v)
            if k in CELL:
                new_cell[k].append(v)
            else:
                raise KeyError('The key {0} you input does not belong to any range!'.format(k))
        return new_control, new_system, new_electrons, new_ions, new_cell

    def write_scf(self, out_file):
        """
        Generate a PWscf input file for given parameters.

        :param out_file: file that you want to write to as output
        :return:
        """
        with open(out_file, 'w') as f:
            f.write("&CONTROL\n")
            f.write("calculation = '{0}'\n".format(self.calc_type))
            for k, v in self.control:
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&SYSTEM\n")
            for k, v in self.system:
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&ELECTRONS\n")
            for k, v in self.eletrons:
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&IONS\n")
            for k, v in self.ions:
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")
            f.write("&CELL\n")
            for k, v in self.cell:
                f.write("{0} = {1}\n".format(k, v))
            f.write("/\n")


class SCFInteractiveGenerator(InteractiveConsole):

