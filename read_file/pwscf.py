#!/usr/bin/env python3
# created at Oct 20, 2017 6:15 PM by Qi Zhang

import collections

from read_file.read_basic import *


class ReadPWscfOutput:
    @staticmethod
    def read_total_energy(inp: str):
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            match = re.findall(
                "!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", f.read())
        return str_list_to_float_list(match)

    KMesh = NamedTuple(
        'KMesh', [('k_grid', List[float]), ('k_shift', List[float])])

    @staticmethod
    def read_k_mesh(inp: str) -> KMesh:
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            for line in f:
                if re.search('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    k_grid = str_list_to_float_list(sp[0:3])
                    k_shift = str_list_to_float_list(sp[3:7])
                    k_mesh = collections.namedtuple(
                        'k_mesh', ['k_grid', 'k_shift'])
                    return k_mesh(k_grid, k_shift)
                else:
                    raise ValueError(
                        "'K_POINTS' not found in your submit! Please check!")
