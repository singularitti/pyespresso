#!/usr/bin/env python3
# created at Aug 18, 2017 11:21 PM by Nil-Zil

from .read_file import ReadOutput


class Compute:
    def __init__(self):
        self.ro = ReadOutput()

    def c_over_a(self, filename: str):
        covera = []
        cp = self.ro.read_final_cell(filename)
        for cell_param in cp:
            a = cell_param[0][0]
            c = cell_param[2][2]
            covera.append(c / a)
        return covera
