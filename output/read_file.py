#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 15:00 by Nil-Zil

from itertools import islice


class ReadOutput:
    def __init__(self, filename: str):
        self.filename = filename

    def read_pv(self) -> tuple:
        """
        Read pressure and volume from file.
        :return: (list, list)
        """
        with open(self.filename, 'r') as f:
            p = []
            v = []
            for line in islice(f, 2, None):
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    break
                else:
                    sp = line.split()
                    p.append(float(sp[0]))
                    v.append(float(sp[1]))
        return p, v

    def read_eos_param(self) -> tuple:
        """
        Read equation of states parameters (volume, bulk modulus and its derivative) from file.
        :return: (float, float, float)
        """
        with open(self.filename, 'r') as f:
            for line in islice(f, 2, None):
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    sp = f.readline().split()
                    v0 = float(sp[2])
                    k0 = float(sp[5])
                    k0p = float(sp[8])
        return v0, k0, k0p

    def read_iter_num(self):
        