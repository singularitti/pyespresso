#!/usr/bin/env python3
# created at Oct 20, 2017 12:49 AM by Qi Zhang

import matplotlib.pyplot as plt

from read_file.read_file import *


class PlotElasticityOutput:
    def __init__(self):
        self.reo = ReadElasticityOutput()

    def plot_cij_vs_pressures(self, inp: str):
        data = self.reo.read_cij_vs_pressures(inp)
        fig, ax = plt.subplots()
        for i, p in enumerate(zip(data.keys(), data.values())):
            for pp in zip([float(p[0])] * len(p[1]), p[1]):
                print(*pp)
                ax.scatter(*pp, label="value" + str(i))
        return fig, ax
