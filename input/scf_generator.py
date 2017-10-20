#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from typing import *
from input.pw_paramlist import *


class SCFGenerator:
    def __init__(self, output_name: str, calc_type: str):
        self.scf_output = output_name
        self.calc_type = calc_type

    @staticmethod
    def readfile(file) -> Dict:
        """
        This method reads input from a given template file. Each line in the file is of form:
            key: value
        or
            key: [value1, value2, ...]
        Then the method generates a dictionary that contains each different keys and corresponding
        value(s).

        :param file: the template file to be read from
        :return: A dictionary that contains each pair of key and value(s) on each line in the file.
        """
        param_dict = {}

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
                k, v = stripped_line.split(':', maxsplit=1)
                k: str = k.strip()
                v: str = v.strip()
                if not k:  # If k is ''
                    raise ValueError('Key is empty string!')
                if not v:  # If v is ''
                    raise ValueError('Value is empty string!')

                if v.startswith('[') and v.endswith(']'):  # v is a list
                    vs = v[1:-2].split(',')  # Skip '[' and ']', then split by ','
                    vx = [x.strip() for x in vs]
                else:
                    vx = v
                param_dict.update({k: vx})

        return param_dict

    def construct_CONTROL(self):
        with open(self.scf_output, 'w') as f:
            f.write("&CONTROL\n")
            f.write("calculation = '%s'\n".format())
