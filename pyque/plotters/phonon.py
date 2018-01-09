#!/usr/bin/env python3
# created at Oct 20, 2017 12:47 AM by Qi Zhang

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from pyque.calculators.phonon import *
from pyque.parsers.elasticity import *


class PhononPlotData:
    def __init__(self, q_path: str):
        """

        :param q_path: A `q_path` is something like 'Γ->M->K->Γ->A->K', indicating the q-point path you are going
            through. Spaces are allowed in this string.
        """
        self.q_path = q_path.upper().replace(' ', '').split('->')
        self.rph = ReadPHononOutput()
        self.cph = ComputePHonon()

    def plot_gnuplot(self, filename: Optional[str] = 'gnuplot'):
        coords_list, bands_list = self.rph.read_gunplot(filename)
        fig, ax = plt.subplots()
        # for i, p in enumerate(zip(coords_list, bands_list)):
        #     ax.plotters(*p, label="band " + str(i))
        ax.plot(coords_list, bands_list)
        # Next we will make some adjustment to plotters all bands.
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
                         box.width, box.height * 0.8])
        fontp = FontProperties()
        fontp.set_size('small')
        ax.legend(loc='center', bbox_to_anchor=(0.5, -0.25), ncol=3, prop=fontp)
        ax.set_xlabel('k-points', fontsize=12)
        ax.set_ylabel("frequency (cm$^{-1})$", fontsize=12)
        ax.set_title('phonon dispersion relation', fontsize=16)

    def generate_plotting_data(self, filename: str, density: IntArray):
        """
        I write this method because I want to separate data generation and plotters steps.

        :param filename:
        :param density:
        :return:
        """
        path_num = len(density)
        qs, bands = self.rph.read_dispersion_relation(filename, density)
        lens = self.cph.q_path_len_list(path_num, qs)
        return qs, bands, lens

    def plot_dos(self, filename: Optional[str] = 'matdyn.dos.out', freq_unit: Optional[str] = 'ev',
                 mode: Optional[str] = 'preview') -> None:
        """
        This method plots the density of states (DOS) of a phonon dispersion relation.

        :param filename: file that contains DOS-data. The val value is 'matdyn.dos.out'.
        :param freq_unit: You can choose from 'ev' or 'cm-1'. The val value is 'ev'.
        :param mode: You can choose to preview/show the plotters ('preview') or save it to file ('save'). The val value
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
        elif freq_unit == 'thz':
            frequency_list = list(map(lambda x: x * 0.02998, frequency_list))
            ax.set_xlabel("frequency (THz)", fontsize=12)
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
            raise ValueError('Unknown readers mode!')


class PlotPhononDispersionRelation(PhononPlotData):
    def dispersion_relation(self, data, color: Optional[str] = 'g', option: Optional[str] = 'hz'):
        """
        This method is used to plotters a phonon dispersion relation for one readers given by matdyn.x.

        :param data: This data is generated by `generate_plotting_data`.
        :param color: Color for the manifold. The val value is green.
        :param option: Specify energy unit, you can submitters 'hz' for hertz, 'thz' for tera-hertz, 'ev' for electron-volt.
            The option is case-insensitive.
        :return:
        """
        qs, bands, lens = data
        fig, axes = plt.subplots(1, len(lens), gridspec_kw={'width_ratios': lens}, sharey='all')

        if re.match('hz', option, flags=re.IGNORECASE):
            option = 'Hz'
            bands = self.cph.frequency_to_hertz(bands)
        elif re.match('thz', option, flags=re.IGNORECASE):
            option = 'THz'
            bands = self.cph.frequency_to_hertz(bands) / 1e12
        elif re.match('ev', option, flags=re.IGNORECASE):
            bands = self.cph.frequency_to_ev(bands)
        else:
            raise ValueError('Unknown option ' + str(option) + ' given!')

        # Main plotters code block
        for i, ax in enumerate(axes):
            x_coordinates = range(len(qs[i]))  # Remap q-points coordinates to x coordinates for plotters
            ax.yaxis.set_ticks_position('none')  # Remove side effect
            ax.plot(x_coordinates, bands[i], color=color)
            ax.get_xaxis().set_ticks([])  # Cancel x-axis ticks
            ax.set_xticks([0, len(qs[i])])  # Specify new ticks position
            ax.set_xticklabels(self.q_path[i])  # Put ticks on x-axis
            ax.set_xlim((min(x_coordinates), max(x_coordinates)))  # To make plotters without inner paddings

        axes[0].set_ylabel(option, fontsize=12)  # Vertical label
        axes[0].yaxis.tick_left()  # Put label on the left

        axes[-1].set_xticks([0, len(qs[-1])])
        axes[-1].set_xticklabels([self.q_path[-2], self.q_path[-1]])
        axes[-1].yaxis.tick_right()
        axes[-1].yaxis.set_ticks_position('right')

        fig.text(0.5, 0.02, 'q-path', fontsize=12)  # Horizontal label
        fig.subplots_adjust(wspace=0, hspace=0)  # Remove spaces between subplots
        return fig, axes

    def plot_multiple_phonon_dispersion(self, files: List[str], density, option):
        colormap = plt.get_cmap('viridis')
        colors: np.ndarray = colormap(np.linspace(0, 1, len(files)))
        qs, bands, lens = self.generate_plotting_data(files[0], density)
        fig, axes = plt.subplots(1, len(lens), gridspec_kw={'width_ratios': lens}, sharey='all')

        font_p = FontProperties()
        font_p.set_size('small')

        for j, file in enumerate(files):
            qs, bands, lens = self.generate_plotting_data(file, density)

            if re.match('hz', option, flags=re.IGNORECASE):
                option = 'Hz'
                bands = self.cph.frequency_to_hertz(bands)
            elif re.match('thz', option, flags=re.IGNORECASE):
                option = 'THz'
                bands = self.cph.frequency_to_hertz(bands) / 1e12
            elif re.match('ev', option, flags=re.IGNORECASE):
                bands = self.cph.frequency_to_ev(bands)
            else:
                raise ValueError('Unknown option ' + str(option) + ' given!')

            for i, ax in enumerate(axes):
                x_coordinates = range(len(qs[i]))  # Remap q-points coordinates to x coordinates for plotters
                ax.yaxis.set_ticks_position('none')  # Remove side effect
                ax.plot(x_coordinates, bands[i], color=colors[j], label=file)
                ax.get_xaxis().set_ticks([])  # Cancel x-axis ticks
                ax.set_xticks([0, len(qs[i])])  # Specify new ticks position
                ax.set_xticklabels(self.q_path[i])  # Put ticks on x-axis
                ax.set_xlim((min(x_coordinates), max(x_coordinates)))  # To make plotters without inner paddings
                box = ax.get_position()
                ax.set_position([box.x0, 0.4, box.width, 0.4])

            # Combine multiple line labels into one legend
            handles, labels = axes[-1].get_legend_handles_labels()
            labels, ids = np.unique(labels, return_index=True)
            handles = [handles[i] for i in ids]
            plt.legend(handles, labels, ncol=1, prop=font_p)

        axes[0].set_ylabel(option, fontsize=12)  # Vertical label
        axes[0].yaxis.tick_left()  # Put label on the left

        axes[-1].set_xticks([0, len(qs[-1])])
        axes[-1].set_xticklabels([self.q_path[-2], self.q_path[-1]])
        axes[-1].yaxis.tick_right()
        axes[-1].yaxis.set_ticks_position('right')

        fig.text(0.5, 0.02, 'q-path', fontsize=12)  # Horizontal label
        fig.subplots_adjust(wspace=0, hspace=0)  # Remove spaces between subplots
        return fig, axes
