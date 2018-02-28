#!/usr/bin/env python3
# created at Jul 19, 2017 16:11 by Qi Zhang

from pyespresso.plotting.plot_basic import *
from pyespresso.lexer.elasticity import *


class PlotVCRelaxOutput(SingleAxes):
    def __init__(self):
        super().__init__()
        self.ro = ReadVCRelaxOutput()
        self.cmap = self.choose_colormap('nipy_spectral')

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
