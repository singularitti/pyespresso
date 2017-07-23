#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 16:35 by Nil-Zil


class EOS(object):
    def __init__(self):
        pass

    @staticmethod
    def vinet(v0: float, k0: float, k0p: float):
        """
        Vinet EOS takes 3 parameters.
        :param v0: float
        :param k0: float
        :param k0p: float
        :return: function
        """

        def p(v: float) -> float:
            """
            Vinet EOS is a function p of v.
            :param v: volume
            :return: float
            """
            r = (v / v0) ** (1 / 3)
            return 3 * k0 * r ** (-2) * (1 - r) * np.exp(3 / 2 * (k0p - 1) * (1 - r))

        return p
