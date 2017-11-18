#!/usr/bin/env python3
# created at Oct 20, 2017 6:17 PM by Qi Zhang

import shlex
from itertools import islice
from operator import itemgetter
from readers.pwscf import *

import numpy as np

from readers.simple_reader import *


class VCRelaxInputfileReader(PWscfInputReader):
    def read_ions_card(self) -> Dict[str, str]:
        return self._read_card('IONS')

    def read_cell_card(self) -> Dict[str, str]:
        return self._read_card('CELL')

    def build_vc_relax_input_tree(self):
        return {'CONTROL': self.read_control_card(),
                'SYSTEM': self.read_system_card(),
                'ELECTRONS': self.read_electrons_card(),
                'IONS': self.read_ions_card(),
                'CELL': self.read_cell_card(),
                'CELL_PARAMETERS': self.read_cell_parameters(),
                'K_POINTS': self.read_k_mesh()}


class VCRelaxOutfileReader:
    @staticmethod
    def read_pv(inp: str) -> Tuple[List[float], List[float]]:
        """
        Read pressure and volume from one file each time.

        :param inp: A single file that is to be read.
        :return: pressures and volumes
        """
        start = 1
        p_list = []
        v_list = []
        with open(inp, 'r') as f:
            for line in islice(f, start, None):
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    break
                else:
                    sp = line.split()
                    p_list.append(float(sp[0]))
                    v_list.append(float(sp[1]))
        return p_list, v_list

    @staticmethod
    def read_eos_param(inp: str) -> Tuple[float, float, float]:
        """
        Read equation of states parameters (volume, bulk modulus and its derivative) from one file each time.

        :param inp: A single file that is to be read.
        :return: zero-pressure volume, zero-pressure bulk modulus, and its first derivative W.R.T. pressure
        """
        v0 = None
        k0 = None
        k0p = None
        start = 1
        count = 2
        with open(inp, 'r') as f:
            for line in islice(f, start, None):
                count += 1
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    sp = f.readline().split()  # Read next line
                    if 'V0' in sp:
                        v0 = float(sp[2])
                    else:
                        print('V0 not found in your file: ' +
                              inp + 'in line: ' + str(count))
                    if 'K0' in sp:
                        k0 = float(sp[5])
                    else:
                        print('K0 not found in your file:' +
                              inp + 'in line: ' + str(count))
                    if 'Kp' in sp:
                        k0p = float(sp[8])
                    else:
                        print("K0' not found in your file:" +
                              inp + 'in line: ' + str(count))
        return v0, k0, k0p

    def _read_eos_from_multiple_files(self, file_list: list, i: int) -> List[float]:
        """
        This method defines the basic logic of the following 3 methods.
        This should not be called directly by users, but wrapped in other methods.

        :param file_list: a list of files
        :param i: index option, {0: V0, 1: K0, 2: K0'}
        :return: a certain EOS parameter
        """
        param_list = []
        for file in file_list:
            param_list.append(self.read_eos_param(file)[i])
        return param_list

    def read_v0_from_files(self, file_list: List[str]) -> List[float]:
        """
        This method reads EOS parameters V0 from each file in a list of files, and collect them as a list.

        :param file_list: a list of files
        :return: a list of volumes read from those files
        """
        return self._read_eos_from_multiple_files(file_list, 0)

    def read_k0_from_files(self, file_list: list) -> List[float]:
        """
        This method reads EOS parameters K0 from each file in a list of files, and collect them as a list.

        :param file_list: a list of files
        :return: a list of bulk modulus read from those files
        """
        return self._read_eos_from_multiple_files(file_list, 1)

    def read_k0p_from_files(self, file_list: list) -> List[float]:
        """
        This method reads EOS parameters K0' from each file in a list of files, and collect them as a list.

        :param file_list: a list of files
        :return: a list of bulk modulus' first derivatives read from those files
        """
        return self._read_eos_from_multiple_files(file_list, 2)

    @staticmethod
    def read_final_cell(inp: str) -> Tuple[List[float], List[np.ndarray]]:
        """
        This method reads result from 'finalcell' file, and collect the data, prepare them for plotters.

        :param inp: A single file that is to be read.
        :return: ([float], [float])
        """
        p_list = []
        cell_params_list = []
        with open(inp, 'r') as f:
            for line in f:
                # If a line starts with '#', it will be regarded as a comment,
                # we do not parse this line.
                fields = shlex.split(line, comments=True)
                if not fields:
                    continue
                if 'Current folder is' in line:
                    p_list.append(
                        float(re.findall("-?[0-9]+\.[0-9]+", line)[0]))
                if 'CELL_PARAMETERS' in line:
                    cell_params = np.zeros((3, 3))
                    for i in range(3):
                        sp = f.readline().split()
                        cell_params[i] = strs_to_floats(sp)
                    cell_params_list.append(cell_params)
        return p_list, cell_params_list

    @staticmethod
    def read_iter_num(inp: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Read iteration number of each test consisting of a series of pressures from one file each time.
        This works if your result is given by checkiternum.sh in this package.

        :param inp: A single file that is to be read.
        :return:
        """
        p_list = []
        iternum_list = []
        with open(inp, 'r') as f:
            for line in islice(f, 0, None):
                p_list.append(float(re.findall("-?\d.*\.\d", line)[0]))
                iternum_list.append(float(f.readline()))
        # Group pressures and iteration numbers first,
        # then sort the latter according to the order of pressures, then
        # un-group new pressures and iteration numbers.
        [p_list, iternum_list] = np.transpose(
            sorted(np.transpose([p_list, iternum_list]), key=itemgetter(0)))
        return p_list, iternum_list
