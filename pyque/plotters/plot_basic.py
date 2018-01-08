#!/usr/bin/env python3
# created at Oct 20, 2017 1:35 AM by Qi Zhang

import matplotlib.pyplot as plt
from typing import *


class SingleAxes:
    def __init__(self):
        self.colormap = self.choose_colormap('viridis')

    @staticmethod
    def choose_colormap(cmap: Optional[str] = 'viridis'):
        return plt.set_cmap(cmap)

    @staticmethod
    def set_legend(lines: Tuple, labels: Tuple):
        plt.legend(lines, labels, loc="best")


class MultipleAxes:
    def __init__(self, nrows, ncols, lens):
        self.fig, self.axes = plt.subplots(nrows, ncols, gridspec_kw={'width_ratios': lens}, sharey='all')
