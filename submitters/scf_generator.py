#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Qi Zhang

from collections import defaultdict
from lazy_property import LazyProperty

from meta.text import TextStream
from typing import DefaultDict
from default_configuerations.generator_mappings import pwscf_generator_config


class PWscfStandardInputCreator(TextStream):
    """

    """

    def read_fixed_input(self) -> DefaultDict[str, str]:
        """
        This method reads submitters from a given template file. Each line in the file is of form:
            key: value
        Then the method generates a dictionary that contains each different keys and corresponding
        value.

        :return: A dictionary that contains each pair of key and value on each line in the file.
        """
        raw_parameters = defaultdict(str)

        for line in self.stream_generator():
            print(line)
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
            if not key:  # if k is '':
                raise ValueError('The key of {0} is an empty string!'.format(value))
            if not value:  # if v is '':
                raise ValueError('The value of {0} is an empty string!'.format(key))
            raw_parameters.update({key: value})

        return raw_parameters

    @LazyProperty
    def raw_parameters(self):
        return self.read_fixed_input()

    def build_parameters(self):
        """
        Split fixed parameters into 5 dictionaries.
        """
        parameters = set()
        for k, v in self.raw_parameters.items():
            parameters.add(pwscf_generator_config(k, v))
        return parameters

    def read_variable_input(self):
        pass
