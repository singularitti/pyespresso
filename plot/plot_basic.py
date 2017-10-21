#!/usr/bin/env python3
# created at Oct 20, 2017 1:35 AM by Qi Zhang

import matplotlib.pyplot as plt
from typing import *


class SingleAxes:
    def __init__(self):
        self.fig, self.ax = plt.subplots()
        self.colormap = self.choose_colormap('viridis')

    @staticmethod
    def choose_colormap(cmap: Optional[str] = 'viridis'):
        return plt.get_cmap(cmap)


class MultipleAxes:
    def __init__(self, nrows, ncols, lens):
        self.fig, self.axes = plt.subplots(nrows, ncols, gridspec_kw={'width_ratios': lens}, sharey='all')
