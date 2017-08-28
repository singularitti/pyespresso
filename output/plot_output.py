#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Nil-Zil

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from output.compute import *


class PlotVCRelax:
    def __init__(self):
        self.ro = ReadVCRelaxOutput()
        self.cmap = plt.get_cmap('nipy_spectral')

    def parse_filenames(self, filename_list):
        pass

    def plot_p_vs_v(self, files: list):
        colors = self.cmap(np.linspace(0, 1, len(files)))
        labels = [file.split(".")[1] for file in files]
        fig, ax = plt.subplots()
        for i, file in enumerate(files):
            [p, v] = self.ro.read_pv(file)
            ax.plot(p, v, 'o-', label=labels[i], color=colors[i])

    def plot_vinet_eos(self, ax):
        # p is a function takes 1 parameter.
        pass

    def plot_v0_vs_temp(self, file_list: List[str]) -> None:
        v0_list = self.ro.read_v0_from_files(file_list)
        temp_list = []
        fig, ax = plt.subplots()
        ax.plot(temp_list, v0_list)
        ax.set_title("$V_0$ vs $T$", fontsize=16)
        ax.set_xlabel("$T$ (K)", fontsize=12)
        ax.set_ylabel("$V_0$ (au$^3$)", fontsize=12)
        plt.show()

    @staticmethod
    def plot_iternum_vs_p(files: list, ax):
        r = ReadVCRelaxOutput()
        for file in files:
            [p, iternum] = r.read_iter_num(file)
            ax.plot(p, iternum)
        ax.set_title('iteration numbers vs pressures on different tests', fontsize=16)
        ax.set_xlabel('pressures (GPa)', fontsize=12)
        ax.set_ylabel('iteration numbers', fontsize=12)


class PlotPHononOutput:
    def __init__(self):
        self.rpb = ReadPHononOutput()
        self.cph = ComputePHonon()

    def plot_gnuplot(self, filename: Optional[str] = 'gnuplot'):
        coords_list, bands_list = self.rpb.read_gunplot(filename)
        fig, ax = plt.subplots()
        for i, (coords, bands) in enumerate(zip(coords_list, bands_list)):
            ax.plot(coords, bands, label="band " + str(i))
        # Next we will make some adjustment to plot all bands.
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
        ks, bands, ls = self.rpb.read_phonon_dispersion('freq.out', 100)
        fig, axes = plt.subplots(nrows=1, ncols=5, sharey='all', gridspec_kw={'width_ratios': ls})
        plt.subplots_adjust(wspace=0, hspace=0)  # Remove spaces between subplots
        for i in range(len(axes)):
            axes[i].get_xaxis().set_ticks([])  # Cancel x-axis ticks
            axes[i].plot(range(len(ks[i])), bands[i])
            axes[i].set_xlim((min(range(len(ks[i]))), max(range(len(ks[i])))))  # To make plot without inner paddings
            axes[i].yaxis.set_ticks_position('none')  # Remove side effect
        fig.suptitle('phonon dispersion relation', fontsize=16)
        axes[0].set_ylabel("frequency (cm$^{-1})$", fontsize=12)
        axes[0].yaxis.tick_left()
        axes[-1].yaxis.tick_right()
        axes[-1].yaxis.set_ticks_position('right')
        plt.show()

    def plot_dos(self, filename: Optional[str] = 'dos.out', freq_unit: Optional[str] = 'ev',
                 mode: Optional[str] = 'preview') -> None:
        """
        This method plots the density of states (DOS) of a phonon dispersion relation.

        :param filename: file that contains DOS-data. The default value is 'dos.out'.
        :param freq_unit: You can choose from 'ev' or 'cm-1'. The default value is 'ev'.
        :param mode: You can choose to preview/show the plot ('preview') or save it to file ('save'). The default value
            is 'preview'.
        :return: None
        """
        fig, ax = plt.subplots()

        frequency_list, dos_list = self.rpb.read_dos(filename)
        if freq_unit == 'ev':
            frequency_list = self.cph.frequency_to_ev(frequency_list)  # Convert unit from cm^-1 to ev
            ax.set_xlabel("electron-volt", fontsize=12)
        elif freq_unit == 'cm-1':
            ax.set_xlabel("frequency (cm$^{-1})$", fontsize=12)
        else:
            raise ValueError('Unknown frequency unit!')

        ax.plot(frequency_list, dos_list)
        ax.set_ylabel('density of states', fontsize=12)
        ax.set_title('phonon density of states', fontsize=16)

        if mode == 'preview':
            plt.show()
        elif mode == 'save':
            fig.savefig("dos" + str(np.random.random_integers(0, 1000)) + ".png", dpi=400)
        else:
            raise ValueError('Unknown output mode!')
