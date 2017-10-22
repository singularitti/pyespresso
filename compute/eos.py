#!/usr/bin/env python3
# created at Jul 19, 2017 16:35 by Qi Zhang

from typing import *

import numpy as np
from scipy.optimize import fsolve


class EOS:
    def vinet(self, v: float, v0: float, k0: float, k0p: float) -> float:
        """
        Vinet EOS takes 3 parameters to form a function P of V.

        :param v: float
        :param v0: float
        :param k0: float
        :param k0p: float
        :return: float
        """
        r = (v / v0) ** (1 / 3)
        return 3 * k0 * r ** (-2) * (1 - r) * np.exp(3 / 2 * (k0p - 1) * (1 - r))

    def solve_vinet(self, p: float, v0: float, k0: float, k0p: float) -> Tuple[np.ndarray, dict, int, str]:
        """
        Suppose you know a certain pressure,
        and you want to find out what is the corresponding volume under this pressure.
        This function will return a numpy array containing the volume.

        :param p:
        :param v0: float
        :param k0: float
        :param k0p: float
        :return:
        """
        def func(v):
            return self.vinet(v, v0, k0, k0p) - p

        return fsolve(func, np.array([100]))


if __name__ == "__main__":
    eos = EOS()
    print(eos.vinet(100, 137.6852, 283.29, 4.86))
    print(eos.solve_vinet(50, 137.6852, 283.29, 4.86))
