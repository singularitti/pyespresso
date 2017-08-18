#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Nil-Zil

import matplotlib.pyplot as plt
import numpy as np

from .read_file import ReadOutput


class PlotOutput:
    def __init__(self):
        self.ro = ReadOutput()

    def plot_p_vs_v(self, files: list, ax):
        cmap = plt.get_cmap('nipy_spectral')
        colors = cmap(np.linspace(0, 1, len(files)))
        labels = [file.split(".")[1] for file in files]
        for i, file in enumerate(files):
            [p, v] = self.ro.read_pv(file)
            ax.plot(p, v, 'o-', label=labels[i], color=colors[i])

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        pass

    def plot_v0_vs_temp(self, files: list) -> list:
        """
        This function reads EOS parameters V0 from each file in files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        v0list = []
        for i in files:
            v0list.append(self.ro.read_eos_param(i)[0])
        return v0list

    def plot_k0_vs_temp(self, files: list) -> list:
        """
        This function reads EOS parameters K0 from each file in files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        k0list = []
        for i in files:
            k0list.append(self.ro.read_eos_param(i)[1])
        return k0list

    def plot_k0p_vs_temp(self, files: list) -> list:
        """
        This function reads EOS parameters K0' from each file in files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        k0plist = []
        for i in files:
            k0plist.append(self.ro.read_eos_param(i)[2])
        return k0plist

    @staticmethod
    def plot_iternum_vs_p(files: list, ax):
        r = ReadOutput()
        for file in files:
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
