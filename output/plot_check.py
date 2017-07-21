#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 16:11 by Nil-Zil

import numpy as np

from eos import EOS
from read_file import ReadOutput


class PlotCheck(object):
    def __init__(self, outputfile: str):
        self.outputfile = outputfile
        self.readoutput = ReadOutput(self.outputfile)
        self.eos = EOS()
        self.p, self.v = self.readoutput.read_pv()
        self.v0, self.k0, self.k0p = self.readoutput.read_eos_param()

    def plot_p_vs_v(self, ax):
        ax.plot(self.v, self.p, 'o-', label="simulation result")

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        p = self.eos.vinet(self.v0, self.k0, self.k0p)
        v = np.linspace(self.v[0], self.v[-1], 200)
        ax.plot(v, list(map(p, v)), label="Vinet EOS")

    @staticmethod
    def plot_labels(ax):
        ax.set_xlabel("volume $(au^3)$", fontsize=12)
        ax.set_ylabel("pressure (GPa)", fontsize=12)
        ax.legend(loc="best")
