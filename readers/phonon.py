#!/usr/bin/env python3
# created at Oct 20, 2017 6:19 PM by Qi Zhang

import numpy as np

from readers.simple_reader import *
from miscellaneous.phonon_params import INPUTPH_card

# Type aliases
IntArray = Union[int, List[int], np.ndarray]


class PhononInputReader(CardReader):
    def __init__(self, in_file):
        super().__init__(in_file, INPUTPH_card)

    def read_inputph_card(self) -> Dict[str, str]:
        return self._read_card('INPUTPH')

    def read_line_of_input(self):
        with open(self.in_file, 'r') as f:
            for line in f:
                if line.strip().startswith('/'):  # Find the end of 'INPUTPH' card
                    line = f.readline()
                    if len(line.strip().split()) == 3:
                        if self.card['ldisp'] == '.false.' and self.card['qplot'] == '.false.':
                            xq = list(map(float, line.split()))
                            self.__dict__['xq'] = xq
                            return xq
                        else:
                            raise ValueError('xq(1), xq(2), xq(3) are given but `ldisp=.true.` or `qplot=.true.`!')
                    elif len(line.strip().split()) == 1:
                        q_points = []
                        if self.card['qplot'] == '.true.':
                            q_points_num = int(line.strip().split()[0])
                            for _ in range(q_points_num):
                                line = f.readline()
                                q_points.append(list(map(float, line.strip().split())))
                                self.__dict__['q-points'] = q_points
                                return q_points
                        else:
                            raise ValueError("Number of q-points is given but `qplot` is not '.true.'!")

    def build_phono_input_tree(self):
        tree = {'INPUTPH': self.read_inputph_card()}
        if hasattr(self, 'xq'):
            tree.update({'xq': self.xq})
        if hasattr(self, 'q-points'):
            tree.update({'q-points': self.q_points})
        return tree

    def __getattr__(self, item):
        if item == 'card':
            self.__dict__['card'] = self.read_inputph_card()
            return self.card
        elif item == 'tree':
            return self.build_phono_input_tree()
        else:
            raise AttributeError('PhononInputReader instance has no attribute {0}'.format(item))


class PhononOutputReader(SingleFileReader):
    def read_dispersion_relation(self, density: IntArray) -> Tuple[np.ndarray, np.ndarray]:
        """
        This method reads phonon dispersion relation returned by matdyn.x.
        This only works for 3-dimensional q-points grid.

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

    def read_q_points(self) -> Dict[str, List[float]]:
        """
        This is exactly same as `_read_reciprocal_points` method, but for phonons, we here call them q-points rather
        than k-points.

        :return: a dictionary of q-points
        """
        return self._read_reciprocal_points()

    def _read_reciprocal_points(self) -> Dict[str, List[float]]:
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
        This method reads through those names and numbers, and set each _name as a key, each 3 k-coordinates as
        its value, forms a dictionary.

        :return: a dictionary
        """
        return self._read_one_column_as_keys(0, lambda x: list(map(float, x)))


class PlotbandReader(SingleFileReader):
    def read_gunplot(self) -> Tuple[List[float], List[float]]:
        """
        Read in coordinates and energy information, and the collect them as an array.
        The file is given by plotband.x, which is an excutable of Quantum ESPRESSO used to plot the phonon dispersion
        relation given by matdyn.x.

        :return: 2 lists of floats, the first is the coordinates in q-space, the second is the energy.
        """
        coordinates, bands = self.read_two_columns()
        return strs_to_floats(coordinates), strs_to_floats(bands)


class DOSReader(SingleFileReader):
    def read_dos(self) -> Tuple[List[float], List[float]]:
        """
        This method reads frequency and density of states from a density of states file generated by matdyn.x
        automatically.

        :return: 2 lists of floats, the first is the energies, the second is the density of states.
        """
        frequencies, dos = self.read_two_columns()  # Tuple[List[str], List[str]]
        return strs_to_floats(frequencies), strs_to_floats(dos)
