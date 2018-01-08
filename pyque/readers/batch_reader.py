#!/usr/bin/env python3
# created at Jan 3, 2018 7:51 PM by Qi Zhang

import re
from itertools import filterfalse
from typing import *
from pyque.readers.simple import SimpleParser


class SlurmBatchReader(SimpleParser):
    def read_shebang(self) -> str:
        """

        :return: '#!/bin/sh' or something like that.
        """
        for line in self.stream_generator():
            if line.startswith('#!'):
                return line

    def read_directives(self) -> Dict[str, str]:
        """
        {'A': 'mphys', 'N': '1', 'tasks-per-node': '24', 'J': 'job', 'time': '02:00:00'}

        :return:
        """

        directives = {}
        for line in self.stream_generator():
            if 'SBATCH' in line.upper():
                directive_name, directive_value = re.match("#\s*SBATCH\s*--?(\w*-?\w*-?\w*)\s?=?(\w+:?\w*:?\w*)",
                                                           line).groups()
                directives.update({directive_name: directive_value})
        return directives

    def read_modules(self) -> Set[str]:
        """
        {'intel-parallel-studio/2017'}

        :return:
        """
        modules = set()
        for line in self.stream_generator():
            if 'module load' in line:
                module, = re.match("module load\s+(\w.*)", line).groups()
                modules.add(module)
        return modules

    def read_commands(self) -> List[str]:
        """
        Except above, the rest lines are all regarded as shell commands.

        :return:
        """
        predicate = lambda s: any(s.startswith(token) for token in ['#', '\n', 'module load'])
        return [line.rstrip("\n") for line in filterfalse(predicate, self.stream_generator())]
