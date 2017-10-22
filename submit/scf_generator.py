#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from typing import *
from submit.pw_paramlist import *
from collections import defaultdict


class SCFGenerator:
    def __init__(self, file: str, output_name: str, calc_type: str):
        self.scf_output = output_name
        self.calc_type = calc_type
        self.parameters = self.read_pwscf_input(file)

    @staticmethod
    def read_pwscf_input(file) -> Dict[str, List[str]]:
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
        parameters: DefaultDict = defaultdict(list)

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
                    raise ValueError('Key is empty string!')
                if not value:  # If v is ''
                    raise ValueError('Value is empty string!')

                if value.startswith('[') and value.endswith(']'):  # If v is a list
                    tmp = value[1:-2].split(',')  # Skip '[' and ']', then split by ','
                    values = [x.strip() for x in tmp]
                else:  # If v is just a single string
                    values = [value]  # Create a list to keep consistency of values
                parameters.update({key: values})

        return parameters

    def combine_dict(self):
        if default_parameters.keys() <= self.parameters.keys():
            raise KeyError('Your given parameters are too less! The calculation cannot be done!')


    def write_control(self):
        with open(self.scf_output, 'w') as f:
            f.write("&CONTROL\n")
            f.write("calculation = '{0}'\n".format(self.calc_type))
            f.write("prefix = '{0}'\n".format(prefix))
