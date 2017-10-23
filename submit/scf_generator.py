#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from typing import *
from submit.pw_paramlist import *
from collections import defaultdict

# Type aliasing
Param = DefaultDict[str, List[str]]


class SCFGenerator:
    def __init__(self, file: str, calc_type: str):
        self.calc_type = calc_type
        self.parameters = self.read_pwscf_input(file)
        self.control, self.system, self.eletrons, self.ions, self.cell = self.split_parameters()

    @staticmethod
    def read_pwscf_input(file) -> Param:
        """
        This method reads submit from a given template file. Each line in the file is of form:
            key: value
        or
            key: [value1, value2, ...]
        Then the method generates a dictionary that contains each different keys and corresponding
        value(s).

        :param file: the template file to be read from
        :return: A dictionary that contains each pair of key and value(s) on each line in the file.
        """
        parameters = defaultdict(list)

        with open(file, 'r') as f:
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

                if value.startswith('[') and value.endswith(']'):  # If v is a list
                    tmp = value[1:-2].split(',')  # Skip '[' and ']', then split by ','
                    values = [x.strip() for x in tmp]
                else:  # If v is just a single string
                    values = [value]  # Create a list to keep consistency of values
                parameters.update({key: values})

        return parameters

    def split_parameters(self) -> Tuple[Param, Param, Param, Param, Param]:
        # If the keys of `default_parameters` are less than that of `self.parameters`
        if default_parameters.keys() <= self.parameters.keys():
            raise KeyError('Your given parameters less than default! The calculation cannot be done!')
        new_control, new_system, new_electrons, new_ions, new_cell = [defaultdict(list)] * 5
        for k, v in self.parameters:
            if k in CONTROL:
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

    def _split_values(self, dictionary):
        pool = defaultdict(list)
        var = defaultdict(list)
        for k, v in dictionary:
            if len(v) == 1:
                pool[k].append(v)
            elif len(v) > 1:
                var[k].append([{k: i} for i in v])
            else:  # len(v) == 0
                raise ValueError('The value of {0} key is empty!'.format(k))
        'for k, v in '


    def write_control(self, file):
        with open(file, 'w') as f:
            f.write("&CONTROL\n")
            f.write("calculation = '{0}'\n".format(self.calc_type))
            for k, v in self.control:
