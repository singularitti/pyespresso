#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 16:11 by Nil-Zil

import numpy as np

from .eos import EOS
from .read_file import ReadOutput


class PlotCheck:
    def __init__(self, outputlist: list, legend: list):
        self.objlist = [ReadOutput(outputfile) for outputfile in outputlist]
        self.eos = EOS()
        self.plists, self.vlists = self.read_multiple_pv()
        self.v0list, self.k0list, self.k0plist = self.read_multiple_eosparam()
        self.legend = legend

    def read_multiple_pv(self):
        pvlists = [obj.read_pv() for obj in self.objlist]
        plists = [pvlist[0] for pvlist in pvlists]
        vlists = [pvlist[1] for pvlist in pvlists]
        return plists, vlists

    def read_multiple_eosparam(self):
        eosparamlist = [obj.read_eos_param() for obj in self.objlist]
        v0list = [eosparam[0] for eosparam in eosparamlist]
        k0list = [eosparam[1] for eosparam in eosparamlist]
        k0plist = [eosparam[2] for eosparam in eosparamlist]
        return v0list, k0list, k0plist

    def plot_p_vs_v(self, ax):
        for i in range(len(self.objlist)):
            ax.plot(self.vlists[i], self.plists[i], 'o-', label=self.legend[i])

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        for i in range(len(self.objlist)):
            p = self.eos.vinet(self.v0list[i], self.k0list[i], self.k0plist[i])
            v = np.linspace(self.vlists[i][0], self.vlists[i][-1], 200)
            ax.plot(v, list(map(p, v)), label=' '.join(['Vinet EOS', self.legend[i]]))

    @staticmethod
    def plot_labels(ax):
        ax.set_xlabel("volume $(au^3)$", fontsize=12)
        ax.set_ylabel("pressure (GPa)", fontsize=12)
        ax.legend(loc="best")
