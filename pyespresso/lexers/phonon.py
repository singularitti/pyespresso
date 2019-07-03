#!/usr/bin/env python3
# created at Oct 20, 2017 6:19 PM by Qi Zhang

from typing import *

import numpy as np

from pyespresso.lexers.simple import SimpleLexer
from pyespresso.tools.strings import strings_to_floats

# Type aliases
IntArray = Union[int, List[int], np.ndarray]


class INPUTPHNamelistLexer:
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_INPUTPH_NAMELIST)


class PHononInputLexer(SimpleLexer):
    def parse_title(self) -> str:
        """
        Read the first non-empty line as the title of the file.

        :return: The title of the file.
        """
        end_index = None
        for i, line in enumerate(self.file_content):
            if '&' in line:
                end_index = i
        return next(s for s in self.file_content[0:end_index] if s)

    def parse_INPUTPH_namelist(self):
        return INPUTPHNamelistLexer(self.infile).lex_namelist()

    def parse_single_q_point(self) -> Optional[np.ndarray]:
        with open(self.infile, 'r') as f:
            for line in f:
                print(line)
                # If a line does not start with '/' (the end of 'INPUTPH' namelist),
                # read next line immediately. If a '/' is met, execute `else` clause.
                if not line.startswith('/'):
                    break
            else:
                # Read next line, and parse the 3 q point coordinates.
                line = f.readline()
                if len(line.strip().split()) == 3:
                    return np.array(strings_to_floats(line.strip().split()))

    def parse_q_points(self):
        with open(self.infile, 'r') as f:
            for line in f:
                if not line.startswith('/'):
                    break
            else:
                line = f.readline()
                if len(line.strip().split()) == 1:
                    q_points_num = int(line.strip().split()[0])
                    q_points = np.empty([q_points_num, 3])
                    for i in range(q_points_num):
                        q_points[i] = strings_to_floats(line.strip().split())
                    return q_points


class PhononOutputLexer(SimpleLexer):
    def read_dispersion_relation(self, density: IntArray) -> Tuple[np.ndarray, np.ndarray]:
        """
        This method reads phonon dispersion relation returned by matdyn.x.
        This only works for 3-dimensional q-points grid.

        :param density: Number of points on each path.
        :return: q-points array and bands array.
        """
        path_num = len(density)
        q_array = np.concatenate([np.zeros([1, density[i], 3]) for i in range(path_num)])
        q = []  # A list of all q-points
        bands = []  # A list of all bands
        with open(self.infile, 'r') as f:
            headline = f.readline()
            nbnd = int(re.findall("nbnd=\s+(\d+)", headline )[0])  # Number of bands for each q-point
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


class PlotbandLexer(SimpleLexer):
    def read_gunplot(self) -> Tuple[Iterable[float], Iterable[float]]:
        """
        Read in coordinates and energy information, and the collect them as an array.
        The file is given by plotband.x, which is an excutable of Quantum ESPRESSO used to plot the phonon dispersion
        relation given by matdyn.x.

        :return: 2 lists of floats, the first is the coordinates in q-space, the second is the energy.
        """
        coordinates, bands = self.read_two_columns()
        return strings_to_floats(coordinates), strings_to_floats(bands)


class DOSLexer(SimpleLexer):
    def read_dos(self) -> Tuple[Iterable[float], Iterable[float]]:
        """
        This method reads frequency and density of states from a density of states file generated by matdyn.x
        automatically.

        :return: 2 lists of floats, the first is the energies, the second is the density of states.
        """
        frequencies, dos = self.read_two_columns()  # Tuple[List[str], List[str]]
        return strings_to_floats(frequencies), strings_to_floats(dos)
