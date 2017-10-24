#!/usr/bin/env python3
# created at Oct 20, 2017 12:52 AM by Qi Zhang

from readers.elasticity import *


class PlotPWscfOutput:
    def __init__(self):
        self.rpw = ReadPWscfOutput()

    def plot_e_vs_k_num(self, filename: str):
        e_list = self.rpw.read_total_energy(filename)
        pass
