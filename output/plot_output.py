#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Nil-Zil

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from miscellaneous.compute import *
from output.generate_test import *


class PlotPWscfOutput:
    def __init__(self):
        self.rpw = ReadPWscfOutput()

    def plot_e_vs_k_num(self, filename: str):
        e_list = self.rpw.read_total_energy(filename)
        pass


class PlotVCRelaxOutput:
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
    def __init__(self, q_path):
        self.q_path = q_path.upper().replace(' ', '').split('->')
        self.rph = ReadPHononOutput()
        self.cph = ComputePHonon()

    def plot_gnuplot(self, filename: Optional[str] = 'gnuplot'):
        coords_list, bands_list = self.rph.read_gunplot(filename)
        print(np.array(coords_list))
        fig, ax = plt.subplots()
        # for i, p in enumerate(zip(coords_list, bands_list)):
        #     ax.plot(*p, label="band " + str(i))
        ax.plot(coords_list, bands_list)
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

    def generate_data_for_plotting(self, filename, density):
        path_num = len(self.q_path) - 1
        qs, bands = self.rph.read_phonon_dispersion(filename, density)
        ls = self.cph.q_path_len_list(path_num, qs)
        return qs, bands, ls

    def plot_phonon_dispersion(self, data, color, option: Optional[str] = 'hz'):
        qs, bands, ls = data
        print(ls)
        fig, axes = plt.subplots(1, 5, gridspec_kw={'width_ratios': ls})

        if option == 'hz':
            bands = self.cph.frequency_to_hertz(bands)
        elif option == 'ev':
            bands = self.cph.frequency_to_ev(bands)
        else:
            raise ValueError('Unknown option' + str(option) + 'given!')

        for i in range(len(axes)):
            axes[i].plot(range(len(qs[i])), bands[i], color=color)
            axes[i].get_xaxis().set_ticks([])  # Cancel x-axis ticks
            axes[i].set_xticks(range(0, 600, 100))
            axes[i].set_xticklabels(self.q_path[i])
            axes[i].set_xlim((min(range(len(qs[i]))), max(range(len(qs[i])))))  # To make plot without inner paddings
            axes[i].yaxis.set_ticks_position('none')  # Remove side effect

            box = axes[i].get_position()
            axes[i].set_position([box.x0, box.y0 + box.height * 0.2,
                                  box.width, box.height * 0.8])
            fontP = FontProperties()
            fontP.set_size('small')
        handles, labels = plt.gca().get_legend_handles_labels()
        labels, ids = np.unique(labels, return_index=True)
        handles = [handles[i] for i in ids]
        axes[2].legend(handles, labels, loc='center', bbox_to_anchor=(0.5, -0.25), ncol=3, prop=fontP)
        axes[2].set_xlabel('q-path', fontsize=12)
        axes[0].set_ylabel('Hz', fontsize=12)
        axes[0].yaxis.tick_left()
        axes[-1].set_xticks([0, 100])
        axes[-1].set_xticklabels([self.q_path[-2], self.q_path[-1]])
        axes[-1].yaxis.tick_right()
        axes[-1].yaxis.set_ticks_position('right')
        plt.subplots_adjust(wspace=0, hspace=0)  # Remove spaces between subplots

    def plot_multiple_phonon_dispersion(self, file_list, density, fig, axes, title):
        cmap = plt.get_cmap('viridis')
        colors = cmap(np.linspace(0, 1, len(file_list)))
        for i, file in enumerate(file_list):
            self.plot_phonon_dispersion(density, axes, file, colors[i])
        fig.suptitle(title, fontsize=16)
        plt.show()

    def plot_dos(self, filename: Optional[str] = 'matdyn.dos.out', freq_unit: Optional[str] = 'ev',
                 mode: Optional[str] = 'preview') -> None:
        """
        This method plots the density of states (DOS) of a phonon dispersion relation.

        :param filename: file that contains DOS-data. The default value is 'matdyn.dos.out'.
        :param freq_unit: You can choose from 'ev' or 'cm-1'. The default value is 'ev'.
        :param mode: You can choose to preview/show the plot ('preview') or save it to file ('save'). The default value
            is 'preview'.
        :return: None
        """
        fig, ax = plt.subplots()
        frequency_list, dos_list = self.rph.read_dos(filename)
        if freq_unit == 'ev':
            frequency_list = self.cph.frequency_to_ev(frequency_list)  # Convert unit from cm^-1 to ev
            ax.set_xlabel("electron-volt", fontsize=12)
        elif freq_unit == 'cm-1':
            ax.set_xlabel("frequency (cm$^{-1})$", fontsize=12)
        else:
            raise ValueError('Unknown frequency unit!')

        ax.set_xlim((min(frequency_list), max(frequency_list)))
        ax.set_ylim((0, max(dos_list)))
        ax.plot(frequency_list, dos_list)
        ax.set_ylabel('density of states', fontsize=12)
        ax.set_title('phonon density of states', fontsize=16)

        if mode == 'preview':
            plt.show()
        elif mode == 'save':
            plt.savefig("dos.pdf")
        else:
            raise ValueError('Unknown output mode!')
