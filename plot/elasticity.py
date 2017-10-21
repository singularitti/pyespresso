#!/usr/bin/env python3
# created at Oct 20, 2017 12:49 AM by Qi Zhang

from miscellaneous.elasticity import *
from plot.plot_basic import *


class PlotElasticityOutput(SingleAxes):
    def __init__(self, filename: str):
        super().__init__()
        self.reo = ElasticityOutputReader(filename)
        self.filename = filename

    def plot_cij_vs_pressures(self):
        data = self.reo.read_cij_vs_pressures()
        for i, p in enumerate(zip(data.keys(), data.values())):
            for pp in zip([float(p[0])] * len(p[1]), p[1]):
                print(*pp)
                self.ax.scatter(*pp, label="value" + str(i))

    def bulk_modulus_voigt_average_vs_pressure(self):
        es = ElasticityCalculator(self.filename)
        voigt = [es.derive_bulk_modulus_voigt_average(e) for e in es.elastic_tensors]
        self.ax.plot(es.pressures, voigt)
