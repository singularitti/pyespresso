#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Nil-Zil

import re
import shlex
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
        start = 1
        p = []
        v = []
        with open(filename, 'r') as f:
            for line in islice(f, start, None):
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
        v0 = None
        k0 = None
        k0p = None
        start = 1
        count = 2
        with open(filename, 'r') as f:
            for line in islice(f, start, None):
                count += 1
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    sp = f.readline().split()  # Read next line
                    if 'V0' in sp:
                        v0 = float(sp[2])
                    else:
                        print('V0 not found in your file: ' + filename + 'in line: ' + str(count))
                    if 'K0' in sp:
                        k0 = float(sp[5])
                    else:
                        print('K0 not found in your file:' + filename + 'in line: ' + str(count))
                    if 'Kp' in sp:
                        k0p = float(sp[8])
                    else:
                        print("K0' not found in your file:" + filename + 'in line: ' + str(count))
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

    def read_v0_from_files(self, files: list) -> list:
        """
        This function reads EOS parameters V0 from each file in a list of files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        v0list = []
        for file in files:
            v0list.append(self.read_eos_param(file)[0])
        return v0list

    def read_k0_from_files(self, files: list) -> list:
        """
        This function reads EOS parameters K0 from each file in a list of files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        k0list = []
        for file in files:
            k0list.append(self.read_eos_param(file)[1])
        return k0list

    def read_k0p_from_files(self, files: list) -> list:
        """
        This function reads EOS parameters K0' from each file in a list of files, and collect them as a list.
        :param files: list(str)
        :return: list(float)
        """
        k0plist = []
        for file in files:
            k0plist.append(self.read_eos_param(file)[2])
        return k0plist

    @staticmethod
    def read_final_cell(filename: str):
        cell_params = []
        with open(filename, 'r') as f:
            for line in f:
                # If line starts with '#', it will be regarded as comment, we do not parse it.
                fields = shlex.split(line, comments=True)
                if not fields:
                    continue
                if 'CELL_PARAMETERS' in line:
                    a = np.zeros((3, 3))
                    for i in range(3):
                        sp = f.readline().split()
                        a[i] = list(map(float, sp))
                    cell_params.append(a)

        return cell_params


class ReadTestCase:
    @staticmethod
    def read_from_ls(filename):
        filelists = []
        with open(filename, 'r') as f:
            for line in f:
                filelists.append(re.split("\n", line)[0])
        return filelists
