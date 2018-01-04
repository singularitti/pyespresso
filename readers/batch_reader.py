#!/usr/bin/env python3
# created at Jan 3, 2018 7:51 PM by Qi Zhang

import re
from itertools import filterfalse
from typing import *


class SlurmBatchReader:
    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None):
        """


        :param instream: If this is given, `infile` argument will be ignored.
        :param infile: If `instream` is not given, this argument will be used.
        """
        if isinstance(instream, str):
            self.instream: str = instream
            self.infile = None
        elif isinstance(infile, str):
            self.infile: str = infile
            self.instream = None
        else:
            raise TypeError('instream and infile cannot be both None! You must specify one of them!')

    def stream_generator(self) -> Iterator[str]:
        if self.instream:  # If `instream` is given and thus not `None`
            for line in self.instream:
                yield line
        else:  # If `infile` is given and thus not `None`
            with open(self.infile, 'r') as f:
                for line in f:
                    yield line

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
