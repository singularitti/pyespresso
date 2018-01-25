#!/usr/bin/env python3
# created at Jul 19, 2017 16:35 by Qi Zhang

from typing import *

import numpy as np
import scipy.optimize


class VinetEoS:
    def __init__(self, v0: float, k0: float, k0p: float):
        """

        :param v0: zero pressure volume of a system
        :param k0: zero pressure bulk modulus of a system
        :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
        """
        self.v0 = v0
        self.k0 = k0
        self.k0p = k0p

    def f_vs_v(self, v: float) -> float:
        """
        The Helmholtz free energy form of VinetEoS EoS, which takes 3 parameters to form a function P of V.

        :param v: the variable
        :return: Helmholtz free energy as a function of V
        """
        r = (v / self.v0) ** (1 / 3)
        xi = 3 / 2 * (self.k0p - 1)
        return 9 * self.k0 * self.v0 / xi ** 2 * (1 + (xi * (1 - r) - 1) * np.exp(xi * (1 - r)))

    def p_vs_v(self, v: float) -> float:
        """
        The pressure form of VinetEoS EoS, which is the first derivative of Helmholtz free energy
        form with respect to V. It takes 3 parameters to form a function P of V.

        :param v: the variable
        :return: pressure as a function of V
        """
        r = (v / self.v0) ** (1 / 3)
        return 3 * self.k0 * r ** (-2) * (1 - r) * np.exp(3 / 2 * (self.k0p - 1) * (1 - r))

    def solve_v_by_p(self, p: float, v0_guess: float) -> Tuple[np.ndarray, dict, int, str]:
        """
        Suppose you know a certain pressure,
        and you want to find out what is the corresponding volume under this pressure.
        This function will return a numpy array containing the volume.

        :param p:
        :param v0_guess: an initial guess for $V_0$
        :return:
        """
        return scipy.optimize.fsolve(lambda v: self.p_vs_v(v) - p, np.array([v0_guess]))

    def fit_p_vs_v(self, vs: List[float], ps: List[float]) -> Tuple[np.ndarray, ...]:
        return scipy.optimize.curve_fit(f=self.p_vs_v, xdata=vs, ydata=ps)


class BirchMurnaghan3rdEoS:
    """
    The third-order Birch–Murnaghan isothermal equation of state
    """

    def __init__(self, v0: float, k0: float, k0p: float):
        """

       :param v0: zero pressure volume of a system
       :param k0: zero pressure bulk modulus of a system
       :param k0p: first derivative of bulk modulus with respect to pressure, evaluated at zero pressure
       """
        self.v0 = v0
        self.k0 = k0
        self.k0p = k0p

    def e_vs_v(self, v: float, *e0) -> float:
        """


        :param v: volume
        :param e0: internal energy at zero pressure
        :return: internal energy
        """
        r = (self.v0 / v) ** (2 / 3)
        if e0:  # If E_0 is given
            e0, = e0
            return e0 + 9 / 16 * self.v0 * self.k0 * ((r - 1) ** 3 * self.k0p + (r - 1) ** 2 * (6 - 4 * r))
        else:
            return 9 / 16 * self.v0 * self.k0 * ((r - 1) ** 3 * self.k0p + (r - 1) ** 2 * (6 - 4 * r))

    def p_vs_v(self, v: float) -> float:
        """

        :param v: volume
        :return: pressure
        """
        r = (self.v0 / v) ** (2 / 3)
        return 3 / 2 * self.k0 * (r ** (7 / 2) - r ** (5 / 2)) * (1 + 3 / 4 * (self.k0p - 4) * (r - 1))
