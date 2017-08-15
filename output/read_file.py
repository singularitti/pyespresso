#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 19, 2017 15:00 by Nil-Zil

import re
from itertools import islice
from operator import itemgetter
import numpy as np


class ReadOutput:
    # def __init__(self, filename):
    #     self.iternum = filename

    @staticmethod
    def read_pv(filename: str) -> tuple:
        """
        Read pressure and volume from file.
        :param filename: str
        :return: (list, list)
        """
        with open(filename, 'r') as f:
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

    @staticmethod
    def read_eos_param(filename) -> tuple:
        """
        Read equation of states parameters (volume, bulk modulus and its derivative) from file.
        :param filename:
        :return: (float, float, float)
        """
        with open(filename, 'r') as f:
            for line in islice(f, 2, None):
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    sp = f.readline().split()
                    v0 = float(sp[2])
                    k0 = float(sp[5])
                    k0p = float(sp[8])
        return v0, k0, k0p

    def read_iter_num(self):
        p = []
        num = []
        with open(self.iternum, 'r') as f:
            for line in islice(f, 0, None):
                p.append(float(re.findall("\-?\d.*\.\d", line)[0]))
                num.append(float(f.readline()))
        pn = np.transpose([p, num])
        pnn = sorted(pn, key=itemgetter(0))
        p = np.transpose(pnn)[0]
        num = np.transpose(pnn)[1]
        return p, num
