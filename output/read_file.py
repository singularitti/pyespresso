#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Nil-Zil

import re
from itertools import islice
from operator import itemgetter
import numpy as np


class ReadOutput:
    @staticmethod
    def read_pv(filename: str) -> tuple:
        """
        Read pressure and volume from one file each time.
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
    def read_eos_param(filename: str) -> tuple:
        """
        Read equation of states parameters (volume, bulk modulus and its derivative) from one file each time.
        :param filename: str
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

    @staticmethod
    def read_iter_num(filename: str) -> tuple:
        """
        Read iteration number of each test consisting of a series of pressures from one file each time.
        This works if your result is given by checkiternum.sh in this package.
        :param filename: str
        :return: (numpy.ndarray, numpy.ndarray)
        """
        p = []
        num = []
        with open(filename, 'r') as f:
            for line in islice(f, 0, None):
                p.append(float(re.findall("-?\d.*\.\d", line)[0]))
                num.append(float(f.readline()))
        # Group p and num first, then sort num according to the order of p, then un-group new p and num.
        [p, num] = np.transpose(sorted(np.transpose([p, num]), key=itemgetter(0)))
        return p, num


class ReadTestCase:
    @staticmethod
    def read_from_ls(filename):
        filelists = []
        with open(filename, 'r') as f:
            for line in f:
                filelists.append(re.split("\n", line)[0])
        return filelists
