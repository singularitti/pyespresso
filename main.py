#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 24, 2017 11:19 AM by Nil-Zil

import matplotlib.pyplot as plt

from output.generate_test import GenerateTest
from output.plot_check import PlotCheck


class Main(object):
    def __init__(self):
        self.test1 = GenerateTest(
            'Out_file', ['750-gpbe', '750-lda', '750-KS'])

    def different_potentials(self, ax):
        filelist = self.test1.generate_filelist()
        legend1 = self.test1.generate_legend()
        test1plot = PlotCheck(filelist, legend1)
        test1plot.plot_p_vs_v(ax)
        # test1plot.plot_vinet_eos(ax)
        test1plot.plot_labels(ax)


if __name__ == '__main__':
    main = Main()
    fig, ax = plt.subplots()
    main.different_potentials(ax)
    plt.show()
