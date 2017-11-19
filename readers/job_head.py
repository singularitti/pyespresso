#!/usr/bin/env python3
# created at Nov 18, 2017 9:13 PM by Qi Zhang

from readers.simple_reader import *


class JobHeadReader(SingleFileReader):
    def read_scheduler(self) -> str:
        return self._match_one_pattern("scheduler:\s*(\w+)", flags=re.IGNORECASE)

    def read_modules(self) -> List[str]:
        with open(self.in_file, 'r') as f:
            match = re.findall('Necessary modules to be loaded.*', f.read(), flags=re.IGNORECASE | re.MULTILINE)
        return list(map(lambda x: x.strip(), match[0].split(':')[1].strip().split(',')))

    def read_nodes_num(self) -> int:
        return self._match_one_pattern("Number of nodes:\s*(\d+)")

    def read_cores_num(self) -> int:
        return self._match_one_pattern("Number of processors:\s*(\d+)")

    def read_shell_shebang(self):
        shebang = []
        with open(self.in_file, 'r') as f:
            for line in f:
                if line.strip().startswith('#'):
                    shebang.append(line)
        return shebang

    def build_job_head_tree(self):
        return {'scheduler': self.read_scheduler(),
                'modules': self.read_modules(),
                'nodes_num': self.read_nodes_num(),
                'cores_num': self.read_cores_num(),
                'shell_shebang': self.read_shell_shebang()}

    def __getattr__(self, item):
        if item == 'tree':
            self.__dict__['tree'] = self.build_job_head_tree()
            return self.tree
        else:
            raise AttributeError("Object has no attribute {0}!".format(item))
