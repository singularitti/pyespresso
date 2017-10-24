#!/usr/bin/env python3
# created at Oct 20, 2017 6:19 PM by Qi Zhang

import numpy as np

from readers.read_basic import *

# Type aliases
IntArray = Union[int, List[int], np.ndarray]


class PHononOutputReader(SimpleReader):
    def read_gunplot(self) -> tuple:
        """
        Read in coordinates and energy information, and the collect them as an array.

        :param inp: A single file that is to be read.
        :return: ([[float]], [[float]])
        """
        coordinates_list = []
        bands_list = []
        with open(self.in_file, 'r') as f:
            for line in f:  # If the line has data
                line = line.strip()
                if line:
                    coordinates_list.append(float(line.split()[0]))
                    bands_list.append(float(line.split()[1]))
        return coordinates_list, bands_list

    def read_dos(self) -> Tuple[List[float], List[float]]:
        """
        This method reads frequency and density of states from a file generated by matdyn.x automatically.

        :return: ([float], [float])
        """
        frequency_list, dos_list = self.read_two_columns()  # Tuple[List[str], List[str]]
        return str_list_to_float_list(frequency_list), str_list_to_float_list(dos_list)

    def read_phonon_dispersion(self, density: IntArray) -> Tuple[np.ndarray, np.ndarray]:
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
        with open(self.in_file, 'r') as f:
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

    def read_k_points(self) -> Dict[str, List[float]]:
        """
        This is exactly same as `_read_reciprocal_points` method, for band structure, we here call them k-points.

        :return: a dictionary
        """
        return self._read_reciprocal_points()

    def read_q_points(self) -> Dict[str, List[float]]:
        """
        This is exactly same as `_read_reciprocal_points` method, but for phonons, we here call them q-points rather
        than k-points.

        :return: a dictionary
        """
        return self._read_reciprocal_points()


class ReadMultiplePHonon(PHononOutputReader):

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
        return [self.read_phonon_dispersion(density) for file in file_list]
