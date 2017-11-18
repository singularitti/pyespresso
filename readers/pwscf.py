#!/usr/bin/env python3
# created at Oct 20, 2017 6:15 PM by Qi Zhang

import os
from collections import namedtuple

import numpy as np

from miscellaneous.pwscf_params import *
from readers.simple_reader import *

# Type alias
KMesh = NamedTuple('KMesh', [('k_grid', List[float]), ('k_shift', List[float])])


class PWscfInputReader(SingleFileReader):
    """
    This class read the pwscf input file in, and parse it to be a tree.
    """

    @staticmethod
    def _section_with_bounds(file, start_pattern, end_pattern) -> Iterator:
        """
        Search in file for the contents between 2 patterns. Referenced from
        [here](https://stackoverflow.com/questions/11156259/how-to-grep-lines-between-two-patterns-in-a-big-file-with-python).

        :param file: file to be read
        :param start_pattern: the pattern labels where the content is going to start, the line contain this pattern is
            ignored
        :param end_pattern: the pattern labels where the content is to an end
        :return: an iterator that can read the file
        """
        section_flag = False
        for line in file:
            if re.match(start_pattern, line, re.IGNORECASE):
                section_flag = True
                line = file.readline()  # If the line is the `start_pattern` itself, we do not parse this line
            if line.startswith(end_pattern):
                section_flag = False
            if section_flag:
                yield line

    def _read_card(self, card_name) -> Dict[str, str]:
        """
        A generic method to read `CONTROL`, `SYSTEM`, `ELECTRONS`, `IONS`, `CELL` cards.

        :param card_name: the card's name, could be 'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', and 'CELL'
        :return: a dictionary that stores the inputted information of the intended card
        """
        card = {}
        start_pattern = '&' + card_name.upper()

        with open(self.in_file, 'r') as f:
            generator = self._section_with_bounds(f, start_pattern, '/')  # '/' separates each card.
            for line in generator:
                stripped_line = line.strip()
                # Use '=' as the delimiter, split the line into key and value.
                key, value = stripped_line.split('=', maxsplit=1)
                key: str = key.strip()
                value: str = value.strip().rstrip(',')  # Ignore trailing comma
                if key in pw_parameters_tree[card_name]:
                    card.update({key: value})
                else:
                    raise KeyError("{0} is not a valid parameter for '{1}' card!".format(key, card_name))

        return card

    def read_control_card(self) -> Dict[str, str]:
        """
        Read everything that falls within 'CONTROL' card.

        :return: a dictionary that stores the inputted information of 'CONTROL' card
        """
        return self._read_card('CONTROL')

    def read_system_card(self) -> Dict[str, str]:
        """
        Read everything that falls within 'SYSTEM' card.

        :return: a dictionary that stores the inputted information of 'CONTROL' card
        """
        return self._read_card('SYSTEM')

    def read_electrons_card(self) -> Dict[str, str]:
        """
        Read everything that falls within 'ELECTRONS' card.

        :return: a dictionary that stores the inputted information of 'ELECTRONS' card
        """
        return self._read_card('ELECTRONS')

    def read_cell_parameters(self) -> np.ndarray:
        """
        Read 3 lines that follows 'CELL_PARAMETERS' string, so there must be no empty line between 'CELL_PARAMETERS' and
        the real cell parameters.

        :return: a numpy array that stores the cell parameters
        """
        cell_param = np.zeros((3, 3))
        with open(self.in_file, 'r') as f:
            for line in f:
                if 'CELL_PARAMETERS' in line.upper():
                    for i in range(3):
                        sp = f.readline().split()
                        cell_param[i] = strs_to_floats(sp)
        return cell_param

    def read_k_mesh(self) -> KMesh:
        """
        Find 'K_POINTS' line in the file, and read the k-mesh.
        If there is no 'K_POINTS' in file, it will read to end, and raise an error.

        :return: a named tuple defined above
        """
        with open(self.in_file, 'r') as f:
            for line in f:
                if re.match('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    grid = strs_to_ints(sp[0:3])
                    shift = strs_to_ints(sp[3:7])
                    k_mesh = namedtuple('k_mesh', ['grid', 'shift'])
                    return k_mesh(grid, shift)
                else:
                    continue
            # Read to EOF
            raise ValueError("'K_POINTS' not found in your input file! Please check!")

    def build_pwscf_input_tree(self) -> Dict[str, Dict[str, Union[str, float, int, NamedTuple, np.ndarray]]]:
        """
        This method combines everything cards, and others in the input file.

        :return: a dictionary that stores every information of the input file
        """
        return {'CONTROL': self.read_control_card(),
                'SYSTEM': self.read_system_card(),
                'ELECTRONS': self.read_electrons_card(),
                'CELL_PARAMETERS': self.read_cell_parameters(),
                'K_POINTS': self.read_k_mesh()}

    def __call__(self) -> dict:
        """
        A method specifies how the class will behave when being called as a function.

        :return: a tree defined by `build_pwscf_input_tree`
        """
        return self.tree

    def __getattr__(self, item):
        """
        This lazily build a `tree` attribute for the class, other attribute except existing ones will be regarded as
        illegal.

        :param item: the attribute user want to get, the only legal one is `tree`
        :return: If the attribute is `tree`, it will return a tree defined by `build_pwscf_input_tree`.
        """
        if item == 'tree':
            self.__dict__['tree'] = self.build_pwscf_input_tree()
            return self.tree
        else:
            raise AttributeError("'{0}' object has no attribute {1}!".format(object, item))

    def __str__(self) -> str:
        """
        A method specifies how the class will behave when being printed to standard output (REPL).

        :return:
        """
        return "The class has a tree like:\n {0}".format(self.tree)

    __repr__ = __str__


class PWscfOutputReader(SingleFileReader):
    def read_lattice_parameter(self) -> float:
        """
        Do not use it "as-is". This is a scaling number, not a 3x3 matrix containing the Cartesian coordinates.

        :return: a scaling number that multiplies `CELL_PARAMETERS` is the read Cartesian coordinates of the primitive
            vectors of a crystal.
        """
        return self._match_only_once("lattice parameter \(alat\)\s+=\s*(\d+\.\d+)", float)

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

    def read_kinetic_energy_cutoff(self) -> float:
        """
        Read kinetic-energy cutoff, at the beginning of the file.

        :return: kinetic-energy cutoff, $E_\text{cut}$
        """
        return self._match_only_once("kinetic-energy cutoff\s+=\s*(-?\d+\.\d+)", float)

    def read_charge_density_cutoff(self) -> float:
        """
        Read charge density cutoff, at the beginning of the file.

        :return: charge density cutoff, $\rho_\text{cut}$
        """
        return self._match_only_once("charge density cutoff\s+=\s*(-?\d+\.\d+)", float)

    def read_atoms_num_per_cell(self) -> int:
        return self._match_only_once("number of atoms\/cell\s+=\s*(\d+)", int)

    def read_atoms_types_num(self) -> int:
        return self._match_only_once("number of atomic types\s+=\s*(\d+)", int)

    def read_nstep(self) -> int:
        return self._match_only_once("nstep\s+=\s*(\d+)", int)

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
                        raise IndexError(error_message)
                    try:
                        for i in range(1, 3):  # Read the rest 2 lines of matrix
                            line = f.readline()
                            stress[i][:] = list(map(float, re.findall(reg2, line)[0][:]))
                    except IndexError:
                        raise IndexError(error_message)
        return stress
