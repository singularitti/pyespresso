#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Qi Zhang
"""
This module will only deal with output files reading processes.
"""

import collections
import re
import shlex
from itertools import islice
from operator import itemgetter
from typing import *

import numpy as np

# Type aliases
IntArray = Union[int, List[int], np.ndarray]


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
        # We add utf-8 support because we may use special characters.
        with open(inp, 'r', encoding='utf-8') as f:
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
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            match = re.findall(
                "!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", f.read())
        return str_list_to_float_list(match)

    KMesh = NamedTuple(
        'KMesh', [('k_grid', List[float]), ('k_shift', List[float])])

    @staticmethod
    def read_k_mesh(inp: str) -> KMesh:
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            for line in f:
                if re.search('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    k_grid = str_list_to_float_list(sp[0:3])
                    k_shift = str_list_to_float_list(sp[3:7])
                    k_mesh = collections.namedtuple(
                        'k_mesh', ['k_grid', 'k_shift'])
                    return k_mesh(k_grid, k_shift)
                else:
                    raise ValueError(
                        "'K_POINTS' not found in your input! Please check!")


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
        # then sort the latter according to the order of pressures, then
        # un-group new pressures and iteration numbers.
        [p_list, iternum_list] = np.transpose(
            sorted(np.transpose([p_list, iternum_list]), key=itemgetter(0)))
        return p_list, iternum_list


class ReadPHononOutput(SimpleRead):
    @staticmethod
    def read_gunplot(inp: str) -> tuple:
        """
        Read in coordinates and energy information, and the collect them as an array.

        :param inp: A single file that is to be read.
        :return: ([[float]], [[float]])
        """
        coordinates_list = []
        bands_list = []
        with open(inp, 'r') as f:
            for line in f:  # If the line has data
                line = line.strip()
                if line:
                    coordinates_list.append(float(line.split()[0]))
                    bands_list.append(float(line.split()[1]))
        return coordinates_list, bands_list

    def read_dos(self, inp: str) -> Tuple[List[float], List[float]]:
        """
        This method reads frequency and density of states from a file generated by matdyn.x automatically.

        :param inp: A single file that is to be read.
        :return: ([float], [float])
        """
        frequency_list, dos_list = self.read_two_columns(inp)  # Tuple[List[str], List[str]]
        return str_list_to_float_list(frequency_list), str_list_to_float_list(dos_list)

    @staticmethod
    def read_phonon_dispersion(inp: str, density: IntArray) -> Tuple[np.ndarray, np.ndarray]:
        """
        This method reads phonon dispersion relation returned by matdyn.x.
        This only works for 3-dimensional q-points grid.

        :param inp: filename
        :param density: Number of points on each path.
        :return: q-points array and bands array.
        """
        path_num = len(density)
        q_array = np.concatenate(
            [np.zeros([1, density[i], 3]) for i in range(path_num)])
        q = []  # A list of all q-points
        bands = []  # A list of all bands
        with open(inp, 'r') as f:
            headline = f.readline()
            nbnd = int(re.findall("nbnd=\s+(\d+)", headline)[0])  # Number of bands for each q-point
            nq = int(re.findall("nks=\s+(\d+)", headline)[0])
            bands_array = np.concatenate(
                [np.zeros([path_num, density[i], nbnd]) for i in range(path_num)])
            for line in f:
                q.append(list(map(float, line.split())))
                newline = f.readline()
                bands.append(list(map(float, newline.split())))

        for i in range(path_num):
            for j in range(density[i]):
                q_array[i][j][:] = q[i * density[i] + j][:]
                bands_array[i][j][:] = bands[i * density[i] + j][:]

        if bands_array.shape[2] == nbnd:
            if q_array.size / 3 == nq:
                return q_array, bands_array
            else:
                raise ValueError('Number of q-points is incorrect! Check your file!')
        else:
            raise ValueError('Number of bands is incorrect! Check your file!')

    def read_multiple_phonon_dispersion(self, file_list: List[str], density: Optional[IntArray]) -> \
            List[Tuple[np.ndarray, np.ndarray]]:
        """
        This method allows you to read a list of q coordinates and a list of bands
        from each file in the list.

        :param file_list: A list of files to be read from.
            You should inspect the same reciprocal-space path for all the files in `file_list`.
        :param density: Used to specify number of points between each 2 neighbor q-points. It should be a list.
        :return:
        """
        return [self.read_phonon_dispersion(file, density) for file in file_list]


class ReadElasticityOutput(SimpleRead):
    @staticmethod
    def read_elastic_tensor(filename: str) -> Tuple[List[float], List[np.ndarray]]:
        """
        Read c_ij file, which looks like:
            Calculated tensor for each pressure:

            P = -10.0
             322.000    123.400    102.250      0.000      0.000      0.000
             123.400    322.000    102.250      0.000      0.000      0.000
             102.250    102.250    307.600      0.000      0.000      0.000
               0.000      0.000      0.000     72.175      0.000      0.000
               0.000      0.000      0.000      0.000     72.175      0.000
               0.000      0.000      0.000      0.000      0.000     99.300

            P = 0.0
             ...

        :param filename: file that stores c_ij matrices
        :return: A list of pressure and a list of 6 by 6 numpy arrays
        """
        pressure_list = []
        elastic_tensor_list = []
        with open(filename, 'r') as f:
            for line in f:
                if not line.strip():  # Ignore blank lines
                    continue
                if 'Calculated tensor for each pressure' in line:
                    f.readline()
                if 'P' in line:
                    pressure_list.append(float(line.split('=')[1]))
                    elastic_tensor = np.empty([6, 6], dtype=np.float64)
                    for i in range(6):
                        elastic_tensor[i] = list(map(float, f.readline().split()))
                    elastic_tensor_list.append(elastic_tensor)
        return pressure_list, elastic_tensor_list

    def read_cij_vs_pressures(self, filename: str) -> Dict[str, List[float]]:
        """
        Read first column as pressures, then the second to the last column as the compliance of the crystal at
        corresponding pressure.

        :param filename:
        :return:
        """
        return self.read_one_column_as_keys(filename, 0, lambda x: list(map(float, x)))
