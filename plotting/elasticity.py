#!/usr/bin/env python3
# created at Oct 20, 2017 12:49 AM by Qi Zhang

import matplotlib.pyplot as plt

from plotting.plot_basic import *
from miscellaneous.elasticity import *


class PlotElasticityOutput(SingleAxisPlot):
    def __init__(self):
        super().__init__()
        self.reo = ReadElasticityOutput()

    def plot_cij_vs_pressures(self, inp: str):
        data = self.reo.read_cij_vs_pressures(inp)
        for i, p in enumerate(zip(data.keys(), data.values())):
            for pp in zip([float(p[0])] * len(p[1]), p[1]):
                print(*pp)
                self.ax.scatter(*pp, label="value" + str(i))

    def bulk_modulus_voigt_average_vs_pressure(self, filename):
        es = ComputeElasticity(filename)
        es.compliance_tensor_list = [es.create_compliance_tensor(e) for e in es.elastic_tensor_list]
        voigt = [es.bulk_modulus_voigt_average(e) for e in es.elastic_tensor_list]
        self.ax.plot(es.pressure_list, voigt)
