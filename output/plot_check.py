#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 16:11 by Nil-Zil

import numpy as np

from .eos import EOS
from .read_file import ReadOutput


class PlotCheck:
    def __init__(self, outputlist: list, *legend: list):
        self.filelist = outputlist
        self.eos = EOS()

    def plot_p_vs_v(self, ax):
        for i in self.filelist:
            [p, v] = ReadOutput().read_pv(i)
            ax.plot(p, v)

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        for i in range(len(self.objlist)):
            p = self.eos.vinet(self.v0list[i], self.k0list[i], self.k0plist[i])
            v = np.linspace(self.vlists[i][0], self.vlists[i][-1], 200)
            ax.plot(v, list(map(p, v)), label=' '.join(['Vinet EOS', self.legend[i]]))

    def plot_v0_vs_temp(self, filelist):
        v0list = []
        for i in filelist:
            v0, k0, k0p = ReadOutput().read_eos_param(i)
            v0list.append(v0)
        return v0list


    def plot_iternum_vs_p(self, ax):
        ls = []
        for i in self.objlist:
            [p, iternum] = i.read_iter_num()
            ax.plot(p, iternum)
            # ls.append(np.transpose([p, iternum]))
            # for j in range(2):
        ax.set_title('iteration numbers vs pressures on different tests', fontsize=16)
        ax.set_xlabel('pressures', fontsize=12)
        ax.set_ylabel('iteration numbers', fontsize=12)

    @staticmethod
    def plot_labels(ax):
        ax.set_xlabel("volume $(au^3)$", fontsize=12)
        ax.set_ylabel("pressure (GPa)", fontsize=12)
        ax.legend(loc="best")
