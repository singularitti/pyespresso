#!/usr/bin/env python3
# created at Oct 19, 2017 11:40 PM by Qi Zhang

from read_file.vcrelax import *


class ComputeVCRelax:
    def __init__(self):
        self.ro = VCRelaxOutputReader()

    def c_over_a(self, filename: str) -> Tuple[List[float], List[float]]:
        """
        This is only used for hexagonal cell.

        :param filename: str
        :return:
        """
        c_over_a_list = []
        p_list, cp_list = self.ro.read_final_cell(filename)
        for cp in cp_list:
            a = cp[0][0]
            c = cp[2][2]
            c_over_a_list.append(c / a)
        return p_list, c_over_a_list
