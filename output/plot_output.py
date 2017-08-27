#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Nil-Zil

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from .read_file import *


class PlotVCRelax:
    def __init__(self):
        self.ro = ReadVCRelaxOutput()

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
        r = ReadVCRelaxOutput()
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


class PlotPHononOutput:
    def __init__(self):
        self.rpb = ReadPHononOutput()

    def plot_gnuplot(self):
        coords, bands = self.rpb.read_gunplot('gnuplot')
        fig, ax = plt.subplots()
        for i, (coord, band) in enumerate(zip(coords, bands)):
            ax.plot(coord, band, label="band " + str(i))

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
                         box.width, box.height * 0.8])
        fontP = FontProperties()
        fontP.set_size('small')
        ax.legend(loc='center', bbox_to_anchor=(0.5, -0.25), ncol=3, prop=fontP)
        ax.set_xlabel('k-points', fontsize=12)
        ax.set_ylabel("frequency (cm$^{-1})$", fontsize=12)
        ax.set_title('phonon dispersion relation', fontsize=16)
        plt.show()

    def plot_phonon_dispersion(self):
        label = ['GM', 'M', 'K', 'GM', 'A', 'K']
        ks, bands, ls = self.rpb.read_phonon_dispersion('freq.out')
        fig, axes = plt.subplots(nrows=1, ncols=5, sharey=True, gridspec_kw={'width_ratios': ls})
        plt.subplots_adjust(wspace=0, hspace=0)  # Remove spaces between subplots
        axes[0].set_ylabel("frequency (cm$^{-1})$", fontsize=12)
        for i in range(len(axes)):
            axes[i].get_xaxis().set_ticks([])
            axes[i].plot(range(len(ks[i])), bands[i])
            axes[i].set_xlim((min(range(len(ks[i]))), max(range(len(ks[i])))))
            axes[i].set_ylim((0, 450))
            axes[i].yaxis.set_ticks_position('none')  # Remove side effect
        plt.suptitle('phonon dispersion relation', fontsize=16)
        # new_tick_label = ____
        # __ax__.set_xticks(new_tick_label)
        # __ax__.set_xticklabels([__ __])
        axes[0].yaxis.tick_left()
        axes[-1].yaxis.tick_right()
        axes[-1].yaxis.set_ticks_position('right')
        plt.show()

    def plot_dos(self):
        frequency, dos = self.rpb.read_dos('dos.out')
        plt.plot(frequency, dos)
        plt.xlabel("frequency (cm$^{-1})$", fontsize=12)
        plt.ylabel('density of states', fontsize=12)
        plt.title('phonon density of states', fontsize=16)
        plt.show()
