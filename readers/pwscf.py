#!/usr/bin/env python3
# created at Oct 20, 2017 6:15 PM by Qi Zhang

import collections
import os
import numpy as np

from readers.read_basic import *

# Type alias
KMesh = NamedTuple('KMesh', [('k_grid', List[float]), ('k_shift', List[float])])


class PWscfOutputReader(SimpleReader):
    def read_total_energy(self) -> float:
        """
        Read total energy from the output file.
        The default unit is Rydberg.

        :return: total energy
        """
        return self._match_only_once("!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", float)

    def read_cell_volume(self) -> float:
        """
        Read the unit-cell volume, which is given at the beginning of the file.
        The default unit is $\text{bohr}^3$.

        :return: unit-cell volume
        """
        return self._match_only_once("unit-cell volume\s+=\s+(\d+\.\d+)", float)

    def read_pressure(self) -> float:
        """
        Read the pressure, which is at the bottom of the file.
        The default unit is kbar.

        :return: pressure at this volume
        """
        return self._match_only_once("P=\s+(-?\d+\.\d+)", float)

    def read_total_stress(self, unit: Optional[str] = 'atomic') -> np.ndarray:
        """
        Read the total stress, in 2 units, from the bottom of the output file.

        :param unit: There are 2 options, 'atomic' and 'kbar', where 'atomic' is the default one.
        :return: a 3x3 numpy array which contains the stress tensor
        """
        stress_atomic = np.zeros((3, 3))
        stress_kbar = np.zeros((3, 3))
        stress = {'atomic': stress_atomic, 'kbar': stress_kbar}
        with open(self.in_file, 'r') as f:
            for line in f:
                if 'total   stress' in line:
                    for i in range(3):  # Read a 3x3 matrix
                        line = f.readline()
                        stress_atomic[i][:] = list(map(float, line.split()))[0:3]
                        stress_kbar[i][:] = list(map(float, line.split()))[3:6]
        return stress[unit]

    def read_kinetic_stress(self) -> np.ndarray:
        return self._read_stress('kinetic')

    def read_local_stress(self) -> np.ndarray:
        return self._read_stress('local')

    def read_nonlocal_stress(self) -> np.ndarray:
        return self._read_stress('nonloc.')

    def read_hartree_stress(self) -> np.ndarray:
        return self._read_stress('hartree')

    def read_exc_cor_stress(self) -> np.ndarray:
        return self._read_stress('exc-cor')

    def read_corecor_stress(self) -> np.ndarray:
        return self._read_stress('corecor')

    def read_ewald_stress(self) -> np.ndarray:
        return self._read_stress('ewald')

    def read_hubbard_stress(self) -> np.ndarray:
        return self._read_stress('hubbard')

    def read_london_stress(self) -> np.ndarray:
        return self._read_stress('london')

    def read_xdm_stress(self) -> np.ndarray:
        return self._read_stress('XDM')

    def read_dft_nl_stress(self) -> np.ndarray:
        return self._read_stress('dft-nl')

    def read_ts_vdw_stress(self) -> np.ndarray:
        return self._read_stress('TS-vdW')

    def read_k_coordinates(self, out_file: str, coordinate_system: Optional[str] = 'crystal'):
        """
        This method can be used to read how many k-points are involved in a PWscf calculation from a file
        outputted by pw.x. Here regular expression is used. Different version of Quantum ESPRESSO may need different
        version of regular expression. If you find a bug, please contact the author at qz2280@columbia.edu.

        :param out_file: output file, where the data read
        :param coordinate_system: Can be 'Cartesian' at any case (upper, lower, etc.); or 'crystal' at any case.
        :return: None
        """
        regexp = "k\(\s+\d+\) = \((\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\), wk =\s+(\d+\.\d+)"
        cart_flag = "cart. coord. in units 2pi/alat"
        cryst_flag = "cryst. coord."

        if coordinate_system.lower() in {'cartesian', 'cart.', 'cart'}:
            system_type = cart_flag
        elif coordinate_system.lower() in {'crystal', 'cryst.', 'cryst'}:
            system_type = cryst_flag
        else:
            raise ValueError(
                "Unknown coordinate system type! It can be either 'Cartesian' or 'crystal'!")

        with open(self.in_file, 'r') as f, open(out_file, 'w') as g:
            flag = False  # If flag is true, read line and match pattern, if not true, no need to match pattern
            k_count = 0  # Count how many k-points have been read
            k_num = None  # How many k-points in total, given by Quantum ESPRESSO
            for line in f:
                if 'number of k points=' in line:
                    k_num = re.findall("number of k points=\s+(\d+)", line)[0]
                    g.write("{0}\n".format(k_num))
                if line.strip() == system_type:
                    flag = True
                if flag:
                    if k_count >= int(k_num):
                        flag = False

                    matches = re.finditer(regexp, line)
                    for match in matches:
                        # The first group is the k-point's coordinates 3-dimensional coordinates, and second
                        # is its weight in first Brillouin zone.
                        for group in match.groups():
                            g.write(group)
                            g.write('   ')  # A separation between first and second group
                        g.write("\n")

                    k_count += 1
        print('Reading done! File is stored at "{0}"'.format(os.path.abspath(out_file)))

    def _read_stress(self, name) -> np.ndarray:
        """
        Read a stress specified by `name`.

        :param name: specify a name of the stress
        :return: a 3x3 numpy array which contains the stress
        """
        reg1 = "\s+stress\s+\(kbar\)\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        reg2 = "(-?\d+\.\d{2})\s*(-?\d+\.\d{2})\s*(-?\d+\.\d{2})"
        error_message = "Some of the numbers in {0} stress are not distinguished!\n".format(name) + \
                        "They maybe too close so there is no space between them! Try to add a space!"
        stress = np.zeros((3, 3))

        with open(self.in_file, 'r') as f:
            for line in f:
                if name in line and 'stress' in line:  # Read the first line of the matrix
                    try:
                        stress[0][:] = list(map(float, re.findall(name + reg1, line)[0][:]))
                    except IndexError:
                        print(error_message)
                        exit()
                    try:
                        for i in range(1, 3):  # Read the rest 2 lines of matrix
                            line = f.readline()
                            stress[i][:] = list(map(float, re.findall(reg2, line)[0][:]))
                    except IndexError:
                        print(error_message)
                        exit()
        return stress


class PWscfInputReader(SimpleReader):
    def read_k_mesh(self) -> KMesh:
        """

        :return: a named tuple defined above
        """
        with open(self.in_file, 'r') as f:
            for line in f:
                if re.search('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    k_grid = str_list_to_float_list(sp[0:3])
                    k_shift = str_list_to_float_list(sp[3:7])
                    k_mesh = collections.namedtuple(
                        'k_mesh', ['k_grid', 'k_shift'])
                    return k_mesh(k_grid, k_shift)
                else:
                    raise ValueError("'K_POINTS' not found in your submit! Please check!")
