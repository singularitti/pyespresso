#!/usr/bin/env python3
# created at Jul 19, 2017 16:35 by Qi Zhang

from typing import *

import numpy as np
from scipy.optimize import fsolve


class Vinet:
    @staticmethod
    def f_vs_v(v: float, v0: float, k0: float, k0p: float) -> float:
        """
        The Helmholtz free energy form of Vinet EoS, which takes 3 parameters to form a function P of V.

        :param v: the variable
        :param v0: zero pressure volume of a system
        :param k0: zero pressure bulk modulus of a system
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        :return: Helmholtz free energy as a function of V
        """
        r = (v / v0) ** (1 / 3)
        xi = 3 / 2 * (k0p - 1)
        return 9 * k0 * v0 / xi ** 2 * (1 + (xi * (1 - r) - 1) * np.exp(xi * (1 - r)))

    @staticmethod
    def p_vs_v(v: float, v0: float, k0: float, k0p: float) -> float:
        """
        The pressure form of Vinet EoS, which is the first derivative of Helmholtz free energy form with respect to V.
        It takes 3 parameters to form a function P of V.

        :param v: the variable
        :param v0: zero pressure volume of a system
        :param k0: zero pressure bulk modulus of a system
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        :return: pressure as a function of V
        """
        r = (v / v0) ** (1 / 3)
        return 3 * k0 * r ** (-2) * (1 - r) * np.exp(3 / 2 * (k0p - 1) * (1 - r))

    def solve_p_vs_v(self, p: float, v0: float, k0: float, k0p: float, v0_guess: Optional[float] = 100) -> \
            Tuple[np.ndarray, dict, int, str]:
        """
        Suppose you know a certain pressure,
        and you want to find out what is the corresponding volume under this pressure.
        This function will return a numpy array containing the volume.

        :param p:
        :param v0: zero pressure volume of a system
        :param k0: zero pressure bulk modulus of a system
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        :param v0_guess: an initial guess for $V_0$
        :return:
        """

        def func(v):
            return self.p_vs_v(v, v0, k0, k0p) - p

        return fsolve(func, np.array([v0_guess]))


class BirchMurnaghan3rd:
    """
    The third-order Birch–Murnaghan isothermal equation of state
    """

    @staticmethod
    def e_vs_v(v: float, v0: float, k0: float, k0p: float, *e0) -> float:
        """


        :param v: volume
        :param v0: volume at zero pressure
        :param k0: bulk modulus at zero pressure
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        :param e0: internal energy at zero pressure
        :return: internal energy
        """
        r = (v0 / v) ** (2 / 3)
        if e0:  # If E_0 is given
            e0, = e0
            return e0 + 9 / 16 * v0 * k0 * ((r - 1) ** 3 * k0p + (r - 1) ** 2 * (6 - 4 * r))
        else:
            return 9 / 16 * v0 * k0 * ((r - 1) ** 3 * k0p + (r - 1) ** 2 * (6 - 4 * r))

    @staticmethod
    def p_vs_v(v: float, v0: float, k0: float, k0p: float) -> float:
        """

        :param v: volume
        :param v0: volume at zero pressure
        :param k0: bulk modulus at zero pressure
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        :return: pressure
        """
        r = (v0 / v) ** (2 / 3)
        return 3 / 2 * k0 * (r ** (7 / 2) - r ** (5 / 2)) * (1 + 3 / 4 * (k0p - 4) * (r - 1))
