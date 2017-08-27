#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Nil-Zil

import re
import shlex
from itertools import islice
from operator import itemgetter
from typing import *

import numpy as np

import miscellaneous.maths as mm


class SimpleRead:
    @staticmethod
    def read_each_line(filename: str) -> List[str]:
        """
        This method reads each line simply from a file, and discard the '\n' character.

        :param filename: A single file that is to be read.
        :return: A list contains all the lines of the file.
        """
        line_list = []
        with open(filename, 'r') as f:
            for line in f:
                line_list.append(re.split('\n', line)[0])
        return line_list

    @staticmethod
    def read_two_columns(filename: str) -> Tuple[List[str], List[str]]:
        """
        This method reads 2 columns from a file.

        :param filename: A single file that is to be read.
        :return: two columns of the file
        """
        col1_list = []
        col2_list = []
        with open(filename, 'r') as f:
            for line in f:
                sp = line.split()
                col1_list.append(sp[0])
                col2_list.append(sp[1])
        return col1_list, col2_list

    @staticmethod
    def read_one_column_as_keys(filename: str, col_index: int, wrapper: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param filename: A single file that is to be read.
        :param col_index: the index of the column that you want to make it as keys.
        :param wrapper: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        key_list = []
        value_list = []
        with open(filename, 'r') as f:
            for line in f:
                sp = line.split()
                key_list.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                value_list.append(wrapper(sp))
        return dict(zip(key_list, value_list))

    def read_k_points(self, filename: str) -> Dict[str, List[float]]:
        """
        Suppose you have a file like this:
            A	0.0000000000	0.0000000000	0.8396948868
            GM	0.0000000000	0.0000000000	0.0000000000
            H	0.8886483087	1.5391840208	0.8396948868
            H2	0.8886483087	1.5391840208   -0.8396948868
            K	0.8886483087	1.5391840208	0.0000000000
            L	1.3329724631	0.7695920104	0.8396948868
            M	1.3329724631	0.7695920104	0.0000000000
        These are the k-points you want to track through.
        This method reads through those names and numbers, and set each name as a key, each 3 k-coordinates as
        its value, forms a dictionary.

        :param filename: file you want to specify your k-points
        :return: a dictionary
        """
        return self.read_one_column_as_keys(filename, 0, lambda x: list(map(float, x)))


class ReadVCRelaxOutput:
    @staticmethod
    def read_pv(filename: str) -> Tuple[List[float], List[float]]:
        """
        Read pressure and volume from one file each time.

        :param filename: A single file that is to be read.
        :return: pressures and volumes
        """
        start = 1
        p_list = []
        v_list = []
        with open(filename, 'r') as f:
            for line in islice(f, start, None):
                if 'Results' in line:  # Read lines until meet "Results for a Vinet EoS fitting"
                    break
                else:
                    sp = line.split()
                    p_list.append(float(sp[0]))
                    v_list.append(float(sp[1]))
        return p_list, v_list

    @staticmethod
    def read_eos_param(filename: str) -> Tuple[float, float, float]:
        """
        Read equation of states parameters (volume, bulk modulus and its derivative) from one file each time.

        :param filename: A single file that is to be read.
        :return: zero-pressure volume, zero-pressure bulk modulus, and its first derivative W.R.T. pressure
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
    def read_final_cell(filename: str) -> Tuple[List[float], List[np.ndarray]]:
        """
        This method reads result from 'finalcell' file, and collect the data, prepare them for plotting.

        :param filename: A single file that is to be read.
        :return: ([float], [float])
        """
        p_list = []
        cell_params_list = []
        with open(filename, 'r') as f:
            for line in f:
                # If a line starts with '#', it will be regarded as a comment, we do not parse it.
                fields = shlex.split(line, comments=True)
                if not fields:
                    continue
                if 'Current folder is' in line:
                    p_list.append(float(re.findall("-?[0-9]+\.[0-9]+", line)[0]))
                if 'CELL_PARAMETERS' in line:
                    cell_params = np.zeros((3, 3))
                    for i in range(3):
                        sp = f.readline().split()
                        cell_params[i] = list(map(float, sp))
                    cell_params_list.append(cell_params)
        return p_list, cell_params_list

    @staticmethod
    def read_iter_num(filename: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Read iteration number of each test consisting of a series of pressures from one file each time.
        This works if your result is given by checkiternum.sh in this package.

        :param filename: A single file that is to be read.
        :return:
        """
        p_list = []
        iternum_list = []
        with open(filename, 'r') as f:
            for line in islice(f, 0, None):
                p_list.append(float(re.findall("-?\d.*\.\d", line)[0]))
                iternum_list.append(float(f.readline()))
        # Group pressures and iteration numbers first,
        # then sort the latter according to the order of pressures, then un-group new pressures and iteration numbers.
        [p_list, iternum_list] = np.transpose(sorted(np.transpose([p_list, iternum_list]), key=itemgetter(0)))
        return p_list, iternum_list


class ReadPHononOutput(SimpleRead):
    @staticmethod
    def read_gunplot(filename: str) -> tuple:
        """
        Read in coordinates and energy information, and the collect them as an array.

        :param filename: A single file that is to be read.
        :return: ([[float]], [[float]])
        """
        coordinates = []
        bands = []
        coordinate = []
        band = []
        with open(filename, 'r') as f:
            for line in f:
                if not line.strip():  # If the line is empty or blank
                    coordinates.append(coordinate)
                    bands.append(band)
                    # coordinate = []  # Clear the list to accept new data
                    # band = []  # Clear the list to accept new data
                else:  # If the line has data
                    coordinate.append(float(line.split()[0]))
                    band.append(float(line.split()[1]))
        return coordinates, bands

    def read_dos(self, filename: str) -> Tuple[List[float], List[float]]:
        """
        This method reads frequency and density of states from a file generated by matdyn.x automatically.

        :param filename: A single file that is to be read.
        :return: ([float], [float])
        """
        frequency_list, dos_list = self.read_two_columns(filename)  # These 2 lists are filled with strings.
        return list(map(float, frequency_list)), list(map(float, dos_list))

    def read_q_points(self, filename: str) -> Dict[str, List[float]]:
        """
        This is exactly same as read read_k_points method, but for phonons, we here call them q-points rather
        than k-points.

        :param filename: a file that you want to specify your q-points
        :return: a dictionary
        """
        return self.read_k_points(filename)

    @staticmethod
    def read_phonon_dispersion(filename: str):
        """
        This method reads phonon dispersion relation returned by matdyn.x.

        :param filename: str
        :return:
        """
        qs = []
        bands = []
        with open(filename, 'r') as f:
            head = f.readline()
            cols = float(re.findall("nbnd=\s+(\d+)", head)[0])
            rows = float(re.findall("nks=\s+(\d+)", head)[0])
            qp = []
            bandp = []
            i = 1
            for line in f:
                qp.append(list(map(float, line.split())))
                newline = f.readline()
                bandp.append(list(map(float, newline.split())))
                if i % 100 == 0:
                    qs.append(qp)
                    bands.append(bandp)
                    qp = []  # Clear
                    bandp = []  # Clear
                i += 1
        qs = np.array(qs)
        bands = np.array(bands)
        if qs.shape[1] * qs.shape[0] == rows and bands.shape[2] == cols:
            ls = []
            for i in range(len(qs)):
                ls.append(mm.compute_3d_distance(qs[i][0], qs[i][-1]))
            return qs, bands, ls
        else:
            raise ValueError('Number of bands or number of k points does not match header!')
