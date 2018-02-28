#!/usr/bin/env python3

import re
import warnings
from itertools import filterfalse
from typing import *

from pyespresso.lexer import simple


class BatchTemplateLexer(simple.SimpleLexer):
    def parse_scheduler(self) -> str:
        return self._match_one_pattern("Scheduler:\s*(\w+)", flags=re.IGNORECASE)

    def parse_modules(self) -> List[str]:
        match = re.findall("Necessary modules.*:\s*(\w.*)", str(self), flags=re.IGNORECASE | re.MULTILINE)
        return list(map(lambda s: s.strip(), match[0].split(',')))

    def parse_nodes_num(self) -> int:
        return self._match_one_pattern("Number of nodes:\s*(\d+)")

    def parse_cores_num(self) -> int:
        return self._match_one_pattern("Number of processors:\s*(\d+)")

    def parse_shebang(self) -> Optional[str]:
        """

        :return: '#!/bin/sh' or something like that.
        """
        for line in self.text_stream.stream_generator():
            if line.startswith('#!'):
                return line.strip("\n")
        else:
            warnings.warn('No shebang was found for this input!')
            return None

    def parse_directives(self):
        for line in self.text_stream.stream_generator():
            if line.startswith('#!'):
                line = next(line)

    def to_dict(self):
        return {'scheduler': self.parse_scheduler(),
                'modules': self.parse_modules(),
                'nodes_num': self.parse_nodes_num(),
                'cores_num': self.parse_cores_num(),
                'shebang': self.parse_shebang()}


class SlurmSystemBatchLexer(simple.SimpleLexer):
    def parse_shebang(self) -> Optional[str]:
        """

        :return: '#!/bin/sh' or something like that.
        """
        for line in self.text_stream.stream_generator():
            if line.startswith('#!'):
                return line.strip("\n")
        else:
            warnings.warn('No shebang was found for this input!')

    def parse_directives(self) -> Dict[str, str]:
        """
        {'A': 'mphys', 'N': '1', 'tasks-per-node': '24', 'J': 'job', 'time': '02:00:00'}

        :return:
        """

        directives = {}
        for line in self.text_stream.stream_generator():
            if 'SBATCH' in line.upper():
                directive_name, directive_value = re.match("#\s*SBATCH\s*--?(\w*-?\w*-?\w*)\s?=?(\w+:?\w*:?\w*)",
                                                           line).groups()
                directives.update({directive_name: directive_value})
        return directives

    def parse_modules(self) -> Set[str]:
        """
        {'intel-parallel-studio/2017'}

        :return:
        """
        modules = set()
        for line in self.text_stream.stream_generator():
            if 'module load' in line:
                module, = re.match("module load\s+(\w.*)", line).groups()
                modules.add(module)
        return modules

    def parse_commands(self) -> List[str]:
        """
        Except above, the rest lines are all regarded as shell commands.

        :return:
        """

        def predicate(s):
            return any(s.startswith(token) for token in ['#', '\n', 'module load'])

        return [line.rstrip("\n") for line in filterfalse(predicate, self.text_stream.stream_generator())]
