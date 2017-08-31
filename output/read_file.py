#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Nil-Zil
"""
This module will only deal with output files reading processes.
"""

import re
import shlex
from itertools import islice
from operator import itemgetter
from typing import *
import collections

import numpy as np

import miscellaneous.maths as mm


def _str_list_to(inp: List[str], to_type) -> List:
    return list(map(to_type, inp))


def str_list_to_int_list(inp: List[str]) -> List[int]:
    return _str_list_to(inp, int)


def str_list_to_float_list(inp: List[str]) -> List[float]:
    return _str_list_to(inp, float)


class SimpleRead:
    @staticmethod
    def read_each_line(inp: str) -> List[str]:
        """
        This method reads each line simply from a file, and discard the '\n' character.

        :param inp: A single file that is to be read.
        :return: A list contains all the lines of the file.
        """
        line_list = []
        with open(inp, 'r') as f:
            for line in f:
                line_list.append(re.split('\n', line)[0])
        return line_list

    @staticmethod
    def read_two_columns(inp: str) -> Tuple[List[str], List[str]]:
        """
        This method reads 2 columns from a file.

        :param inp: A single file that is to be read.
        :return: two columns of the file
        """
        col1_list = []
        col2_list = []
        with open(inp, 'r') as f:
            for line in f:
                if not line.split():
                    f.readline()
                else:
                    sp = line.split()
                    col1_list.append(sp[0])
                    col2_list.append(sp[1])
        return col1_list, col2_list

    @staticmethod
    def read_one_column_as_keys(inp: str, col_index: int, wrapper: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param inp: A single file that is to be read.
        :param col_index: the index of the column that you want to make it as keys.
        :param wrapper: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        key_list = []
        value_list = []
        with open(inp, 'r', encoding='utf-8') as f:  # We add utf-8 support because we may use special characters.
            for line in f:
                sp = line.split()
                key_list.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                value_list.append(wrapper(sp))
        return dict(zip(key_list, value_list))

    def _read_reciprocal_points(self, inp: str) -> Dict[str, List[float]]:
        """
        Suppose you have a file like this:
            A	0.0000000000	0.0000000000	0.5000000000
            Γ   0.0000000000	0.0000000000	0.0000000000
            H	0.3333333333	0.3333333333	0.5000000000
            H2	0.3333333333	0.3333333333   -0.5000000000
            K	0.3333333333	0.3333333333	0.0000000000
            L	0.5000000000	0.0000000000	0.5000000000
            M	0.5000000000	0.0000000000	0.0000000000
        These are the k-points you want to track through.
        This method reads through those names and numbers, and set each name as a key, each 3 k-coordinates as
        its value, forms a dictionary.

        :param inp: file you want to specify your k-points
        :return: a dictionary
        """
        return self.read_one_column_as_keys(inp, 0, lambda x: list(map(float, x)))

    def read_k_points(self, inp: str) -> Dict[str, List[float]]:
        """
        This is exactly same as read _read_reciprocal_points method, for band structure, we here call them k-points.

        :param inp: a file that you want to specify your q-points
        :return: a dictionary
        """
        return self._read_reciprocal_points(inp)

    def read_q_points(self, inp: str) -> Dict[str, List[float]]:
        """
        This is exactly same as read _read_reciprocal_points method, but for phonons, we here call them q-points rather
        than k-points.

        :param inp: a file that you want to specify your q-points
        :return: a dictionary
        """
        return self._read_reciprocal_points(inp)


class ReadPWscfOutput:
    @staticmethod
    def read_total_energy(inp: str):
        with open(inp, 'r') as f:
            match = re.findall(
                "!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", f.read())
        return str_list_to_float_list(match)

    KMesh = NamedTuple('KMesh', [('k_grid', List[float]), ('k_shift', List[float])])

    @staticmethod
    def read_k_mesh(inp: str) -> KMesh:
        with open(inp, 'r') as f:
            for line in f:
                if re.search('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    k_grid = str_list_to_float_list(sp[0:3])
                    k_shift = str_list_to_float_list(sp[3:7])
                    k_mesh = collections.namedtuple('k_mesh', ['k_grid', 'k_shift'])
                    return k_mesh(k_grid, k_shift)
                else:
                    raise ValueError("'K_POINTS' not found in your input! Please check!")


class ReadVCRelaxOutput:
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
        This method reads result from 'finalcell' file, and collect the data, prepare them for plotting.

        :param inp: A single file that is to be read.
        :return: ([float], [float])
        """
        p_list = []
        cell_params_list = []
        with open(inp, 'r') as f:
            for line in f:
                # If a line starts with '#', it will be regarded as a comment, we do not parse it.
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
                        cell_params[i] = str_list_to_float_list(sp)
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
        # then sort the latter according to the order of pressures, then un-group new pressures and iteration numbers.
        [p_list, iternum_list] = np.transpose(
            sorted(np.transpose([p_list, iternum_list]), key=itemgetter(0)))
        return p_list, iternum_list


class ReadPHononOutput(SimpleRead):
    def read_gunplot(self, inp: str) -> tuple:
        """
        Read in coordinates and energy information, and the collect them as an array.

        :param inp: A single file that is to be read.
        :return: ([[float]], [[float]])
        """
        coordinates_list, bands_list = self.read_two_columns(inp)
        return str_list_to_float_list(coordinates_list), str_list_to_float_list(bands_list)

    def read_dos(self, inp: str) -> Tuple[List[float], List[float]]:
        """
        This method reads frequency and density of states from a file generated by matdyn.x automatically.

        :param inp: A single file that is to be read.
        :return: ([float], [float])
        """
        frequency_list, dos_list = self.read_two_columns(inp)  # Tuple[List[str], List[str]]
        return str_list_to_float_list(frequency_list), str_list_to_float_list(dos_list)

    @staticmethod
    def read_density_on_path(inp: str):
        return [100] * 5

    def read_phonon_dispersion(self, inp: str, q_path):
        """
        This method reads phonon dispersion relation returned by matdyn.x.

        :param inp:
        :param q_path:
        :return:
        """
        dens = self.read_density_on_path('aa')
        qp = q_path.upper().replace(' ', '').split('->')
        num_of_paths = len(qp) - 1
        q_list = np.concatenate(
            [np.zeros([num_of_paths, dens[i], 3]) for i in range(num_of_paths)])
        q = []  # A list of all q-points
        bands = []  # A list of all bands
        with open(inp, 'r') as f:
            headline = f.readline()
            # Number of bands for each q-point
            nbnd = int(re.findall("nbnd=\s+(\d+)", headline)[0])
            nks = int(re.findall("nks=\s+(\d+)", headline)[0])
            bands_list = np.concatenate(
                [np.zeros([num_of_paths, dens[i], nbnd]) for i in range(num_of_paths)])
            for line in f:
                q.append(list(map(float, line.split())))
                newline = f.readline()
                bands.append(list(map(float, newline.split())))

        for i in range(num_of_paths):
            for j in range(dens[i]):
                q_list[i][j][:] = q[i * dens[i] + j][:]
                bands_list[i][j][:] = bands[i * dens[i] + j][:]

        q_path_len_list = []
        for i in range(num_of_paths):
            q_path_len_list.append(
                mm.compute_3d_distance(q_list[i][0], q_list[i][-1]))

        return q_list, bands_list, q_path_len_list

        # Reading file is finished
        # if cols is None or rows is None:
        #     raise NameError("'nbnd' or 'nks' not found in your file! Please check it!")
        # elif q_list.shape[0] * q_list.shape[1] == rows and bands_list.shape[2] == cols:
        #     return q_list, bands_list
        # else:
        #     try:
        #         print(q_list.shape[1] * q_list.shape[0], rows)
        #         print(bands_list.shape[2], cols)
        #     except ValueError:
        #         raise ValueError('Number of bands or number of k points does not match header!')
