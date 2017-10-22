#!/usr/bin/env python3
# created at Aug 18, 2017 11:21 PM by Qi Zhang

import compute.maths as mm
from compute.converter import *
from read_file.elasticity import *

# Type aliases
IntArray = Union[int, List[int], np.ndarray]


class PathGenerator:
    @staticmethod
    def linspace_3d(point1: List[float], point2: List[float], dens: int, **option) -> np.ndarray:
        """
        This method generates a series of evenly-spaced points between 2 given points.

        :param point1: point 1
        :param point2: point 2
        :param dens: density, i.e., number of points between 2 given points. If option is endpoint=True,
            then dens includes the start point and the end point.
        :param option: It is the same as options for np.linspace. I suggest that you use endpoint=True,
            which is the default for np.linspace.
        :return: a series of evenly-spaced points
        """
        dist: np.float64 = mm.compute_3d_distance(point1, point2)
        normalize_vec = (np.array(point2) - np.array(point1)) / dist
        return np.array(point1) + [x * normalize_vec for x in np.linspace(0, dist, dens, **option)]

    @staticmethod
    def check_density(path_num: int, density: IntArray) -> IntArray:
        """
        If you give an `density` argument, it will check what actually `density` is: an integer, a list of integers,
        or an numpy array of integers. I grab this function out because it is constantly reused.

        :param path_num: This specify how many paths do you want to have. Like, in 'Γ->M->K->Γ->A->K', you have 6 points,
            thus 5 paths.
        :param density: Number of points on each path
        :return: If the type of `density` is correct, return `density` back.
        """
        # Check if density is given a reasonable type, the only allowed are int, list, and np.ndarray.
        if isinstance(density, int):
            density = [density] * path_num
        elif isinstance(density, list) or isinstance(density, np.ndarray):
            if len(density) != path_num:  # The length of density should be the same to number of q-points minus 1.
                raise ValueError('The length of k-point density is incorrect!')
        else:
            raise TypeError('The type of k-point density is incorrect!')
        return density


class ReciprocalPathGenerator(PathGenerator):
    def __init__(self, inp: str, reci_path: str):
        """
        If you give a path like 'Γ->M->K->Γ->A->K', then this method will automatically generate---according
        to your file given, which records each q-points' coordinate---the coordinate of each point in your path.
        The file looks like:
            A	0.0000000000	0.0000000000	0.5000000000
            Γ   0.0000000000	0.0000000000	0.0000000000
            H	0.3333333333	0.3333333333	0.5000000000
            H2	0.3333333333	0.3333333333   -0.5000000000
            K	0.3333333333	0.3333333333	0.0000000000
            L	0.5000000000	0.0000000000	0.5000000000
            M	0.5000000000	0.0000000000	0.0000000000
        See read_file.read_file.SimpleRead._read_reciprocal_points for more information.

        :param inp: As explained above.
        :param reci_path:  A specific q-path you are interested in, like 'GM->M->K->GM->A->K'.
            Each q-point should be separated by a '->' character, spaces are allow, other characters are not allowed.
            Special character, like 'Γ' in utf-8 is allowed.
            Note that if you use 'Γ' in your `inp` file, you should also use 'Γ' in `reci_path`!
        """
        self.rph = ReadPHononOutput()
        self.reci_dict: Dict[str, List[float]] = self.rph.read_q_points(inp)
        self.reci_path: List[str] = reci_path.upper().replace(' ', '').split('->')

    def generate_k_path(self, density: Optional[IntArray], out: str) -> np.ndarray:
        """
        This is just a wrapper for `generate_reciprocal_path`, with writing to file function.

        :param density: Used to specify number of points between each 2 neighbor q-points. If it is an integer,
            the method will automatically generate an array filled with same value. If it is already an array or list,
            the length of them should be the number of q-points minus 1. The default value is 100.
        :param out: filename you want to write to. The coordinates of the k-points on the path will be written.
        :return: the coordinates of the k-points on the path.
        """
        coords = self._generate_reciprocal_path(density)
        self._write_path_to_file(coords, out)
        return coords

    def generate_q_path(self, density: Optional[IntArray], out: str) -> np.ndarray:
        """
        This is just a wrapper for `generate_reciprocal_path`, with writing to file function.

        :param density: Used to specify number of points between each 2 neighbor q-points. If it is an integer,
            the method will automatically generate an array filled with same value. If it is already an array or list,
            the length of them should be the number of q-points minus 1. The default value is 100.
        :param out: filename you want to write to. The coordinates of the q-points on the path will be written.
        :return: the coordinates of the q-points on the path.
        """
        coords = self._generate_reciprocal_path(density)
        self._write_path_to_file(coords, out)
        return coords

    def _generate_reciprocal_path(self, density: Optional[IntArray]) -> np.ndarray:
        """
        This method will generate an array of coordinates of points on the path, specified in `self.__init__` function,
        making use of 3D linear interpolation of a straight line.

        :param density: Used to specify number of points between each 2 neighbor q-points. If it is an integer,
            the method will automatically generate an array filled with same value. If it is already an array or list,
            the length of them should be the number of q-points minus 1.
        :return: a numpy array of shape `path_num[i], dens[i], 3`, where i ranges from 0 to `path_num - 1`.
        """
        path_num: int = len(self.reci_path) - 1  # Number of paths

        density = self.check_density(path_num, density)

        reci_coords = []
        for i in range(path_num):
            reci_coords.append(
                self.linspace_3d(self.reci_dict[self.reci_path[i]], self.reci_dict[self.reci_path[i + 1]], density[i]))
        return np.array(reci_coords)

    def _write_path_to_file(self, coords: np.ndarray, out: str) -> None:
        """
        Write the above result to file `out`.

        :param coords: As explained above.
        :param out: read_file filename.
        :return: None
        """
        with open(out, 'wb') as f:
            f.write(('The reciprocal path is: ' + '->'.join(self.reci_path) +
                     ", and the number of points is: \n").encode('utf-8'))
            f.write((str(coords.size / 3) + "\n").encode('utf-8'))
            for row in coords:
                np.savetxt(f, row)


class ComputePHonon:
    def __init__(self):
        self.rpb = ReadPHononOutput()

    @staticmethod
    def frequency_to_ev(frequency_list: Union[int, float, List, np.ndarray]) -> np.ndarray:
        """
        This method converts the frequency read from density of states calculation read_file to electron-volt.

        :return: energy in unit of electron-volt
        """
        return call_simple_converter('e', frequency_list, 'cm-1', 'ev')

    @staticmethod
    def frequency_to_hertz(frequency_list: Union[int, float, List, np.ndarray]) -> np.ndarray:
        """
        This method converts the frequency read from density of states calculation read_file to hertz.
        
        :param frequency_list: 
        :return: energy in unit of hertz
        """
        return call_simple_converter('e', frequency_list, 'cm-1', 'hz')

    @staticmethod
    def q_path_len_list(path_num, q_list):
        q_path_len_list = []
        for i in range(path_num):
            q_path_len_list.append(
                mm.compute_3d_distance(q_list[i][0], q_list[i][-1]))
        return q_path_len_list
