#!/usr/bin/env python3
# created at Oct 20, 2017 1:35 AM by Qi Zhang

import matplotlib.pyplot as plt


class SingleAxisPlot:
    def __init__(self):
        self.fig, self.ax = plt.subplots()

    def set_title(self, title: str):
        """
        Convert matplotlib homonymic API to native API.

        :param title: figure title
        :return: None
        """
        self.ax.set_title(title)

    def legend(self):
        """
        Convert matplotlib homonymic API to native API.

        :return: None
        """
        self.ax.legend(loc="best")

    def show(self):
        """
        Convert matplotlib homonymic API to native API.

        :return: None
        """
        plt.show()

    def savefig(self, path):
        """
        Convert matplotlib homonymic API to native API.

        :return: None
        """
        self.fig.savefig(path)
