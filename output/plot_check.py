#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Nil-Zil

import matplotlib.pyplot as plt
import numpy as np

from .read_file import ReadOutput


class PlotOutput:
    @staticmethod
    def plot_p_vs_v(filelist, ax):
        cmap = plt.get_cmap('nipy_spectral')
        colors = cmap(np.linspace(0, 1, len(filelist)))
        labels = [file.split(".")[1] for file in filelist]
        for i, file in enumerate(filelist):
            [p, v] = ReadOutput().read_pv(file)
            ax.plot(p, v, 'o-', label=labels[i], color=colors[i])

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        for i in range(len(self.objlist)):
            p = self.eos.vinet(self.v0list[i], self.k0list[i], self.k0plist[i])
            v = np.linspace(self.vlists[i][0], self.vlists[i][-1], 200)
            ax.plot(v, list(map(p, v)), label=' '.join(['Vinet EOS', self.legend[i]]))

    @staticmethod
    def plot_v0_vs_temp(filelist):
        v0list = []
        for i in filelist:
            v0, k0, k0p = ReadOutput().read_eos_param(i)
            v0list.append(v0)
        return v0list

    @staticmethod
    def plot_k0_vs_temp(filelist):
        k0list = []
        for i in filelist:
            v0, k0, k0p = ReadOutput().read_eos_param(i)
            k0list.append(k0)
        return k0list

    @staticmethod
    def plot_k0p_vs_temp(filelist):
        k0plist = []
        for i in filelist:
            v0, k0, k0p = ReadOutput().read_eos_param(i)
            k0plist.append(k0p)
        return k0plist

    @staticmethod
    def plot_iternum_vs_p(filelist, ax):
        r = ReadOutput()
        for file in filelist:
            [p, iternum] = r.read_iter_num(file)
            ax.plot(p, iternum)
        ax.set_title('iteration numbers vs pressures on different tests', fontsize=16)
        ax.set_xlabel('pressures (GPa)', fontsize=12)
        ax.set_ylabel('iteration numbers', fontsize=12)

    @staticmethod
    def plot_labels(ax):
        ax.set_xlabel("volume $(au^3)$", fontsize=12)
        ax.set_ylabel("pressure (GPa)", fontsize=12)
        ax.legend(loc="best")
