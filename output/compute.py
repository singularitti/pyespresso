#!/usr/bin/env python3
# created at Aug 18, 2017 11:21 PM by Nil-Zil

from miscellaneous.converter import *
from output.read_file import *


class ComputeVCRelax:
    def __init__(self):
        self.ro = ReadVCRelaxOutput()

    def c_over_a(self, filename: str) -> Tuple[List[float], List[float]]:
        """
        This is only used for hexagonal cell.

        :param filename: str
        :return:
        """
        c_over_a_list = []
        p_list, cp_list = self.ro.read_final_cell(filename)
        for cp in cp_list:
            a = cp[0][0]
            c = cp[2][2]
            c_over_a_list.append(c / a)
        return p_list, c_over_a_list


class ComputePHonon:
    def __init__(self):
        self.rpb = ReadPHononOutput()

    @staticmethod
    def frequency_to_ev(frequency_list: List[float]) -> List[float]:
        """
        This method converts the frequency read from density of states calculation output to electron-volt.

        :return: energy in unit of electron-volt
        """
        return call_simple_converter('e', frequency_list, 'cm-1', 'ev')

    def generate_q_path(self, filename: str, output_filename: str, q_path: str,
                        dens: Union[None, int, List[int], np.ndarray] = 100,
                        mode: Optional[str] = 'default') -> None:
        """
        If you give a path like 'GM->M->K->GM->A->K', then this method will automatically generate---according
        to your file given, which records each q-points' coordinate---the coordinate of each point in your path.
        The file is like:
            A	0.0000000000	0.0000000000	0.8396948868
            Γ   0.0000000000	0.0000000000	0.0000000000
            H	0.8886483087	1.5391840208	0.8396948868
            H2	0.8886483087	1.5391840208   -0.8396948868
            K	0.8886483087	1.5391840208	0.0000000000
            L	1.3329724631	0.7695920104	0.8396948868
            M	1.3329724631	0.7695920104	0.0000000000
        See output.read_file.SimpleRead.read_k_points for more information.

        :param filename: As explained above.
        :param output_filename: Write result to here.
        :param q_path: A specific q-path you are interested in, like 'GM->M->K->GM->A->K'. Each q-point should be
            separated by a '->' character, spaces are allow, other characters are not allowed.
        :param dens: Used to specify number of points between each 2 neighbor q-points. If it is an integer,
            the method will automatically generate an array filled with same value. If it is already an array or list,
            the length of them should be the number of q-points minus 1. The default value is 100.
        :param mode: If you are in debugging mode, then will not write to file but directly print the result,
            if you are in default mode, write result to file, if you input a wrong mode, error will be raised. The
            default value is 'default'.
        :return: None
        """
        qp = q_path.upper().replace(' ', '').split('->')
        path_segment = len(qp) - 1

        # Check if dens is given a reasonable value
        if isinstance(dens, int):
            density = [dens] * path_segment
        elif isinstance(dens, list) or isinstance(dens, np.ndarray):
            if not len(dens) == path_segment:
                raise ValueError('The length of k-point density is incorrect!')
        else:
            raise TypeError('The type of k-point density is incorrect!')

        q_dict = self.rpb.read_q_points(filename)
        coords = []
        for i in range(path_segment):
            coords.append(self.generate_3d_segment(
                q_dict[qp[i]], q_dict[qp[i + 1]], density[i]))
        coords = np.array(coords).reshape((path_segment * dens, 3))

        if mode == 'debug':
            print(coords)
        elif mode == 'default':
            with open(output_filename, 'ab') as f:
                f.write(('K path is: ' + q_path).encode('utf-8'))
                f.write(str(coords.size).encode('utf-8'))
                np.savetxt(output_filename, coords, fmt='%f')
        else:
            raise ValueError('You input a wrong mode!')

    @staticmethod
    def generate_3d_segment(point1: List[float], point2: List[float], dens: int) -> np.ndarray:
        """
        This method generates a series of evenly-spaced points between 2 given points.
        :param point1: point 1
        :param point2: point 2
        :param dens: density, i.e., number of points between 2 given points
        :return: a series of evenly-spaced points
        """
        l = mm.compute_3d_distance(point1, point2)
        normalize_vec = (np.array(point2) - np.array(point1)) / l
        return np.array(point1) + \
               [x * normalize_vec for x in np.linspace(0, l, dens)]
