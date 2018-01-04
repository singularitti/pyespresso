#!/usr/bin/env python3
# created at Oct 20, 2017 6:15 PM by Qi Zhang

from basics.pwscf import *
from readers.simple import *
from basics.pw_params import *

# Type alias
KPoints = NamedTuple('KPoints', [('grid', List[float]), ('offsets', List[float])])

AtomicSpecies = namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])


def write_to_file(obj: object, out_file: str):
    if isinstance(obj, PWStandardInput):
        obj.write_to_file(out_file)
    else:
        raise TypeError('Input object is not an {0}!'.format('SCFStandardInput'))


# ====================================== The followings are input readers. ======================================
class CONTROLNamelistParser(NamelistParserGeneric):
    def __init__(self, infile):
        super().__init__(infile, CONTROL_namelist)


class SYSTEMNamelistParser(NamelistParserGeneric):
    def __init__(self, infile):
        super().__init__(infile, SYSTEM_namelist)


class ELECTRONSNamelistParser(NamelistParserGeneric):
    def __init__(self, infile):
        super().__init__(infile, ELECTRONS_namelist)


class IONSNamelistParser(NamelistParserGeneric):
    def __init__(self, infile):
        super().__init__(infile, IONS_namelist)


class CELLNamelistParser(NamelistParserGeneric):
    def __init__(self, infile):
        super().__init__(infile, CELL_namelist)


class PWInputParser(SingleFileParser):
    """
    This class read an scf input file in, and parse it to be a tree.
    """

    def parse_CONTROL_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'CONTROL' card.

        :return: a dictionary that stores the inputted information of 'CONTROL' card
        """
        return CONTROLNamelistParser(self.infile).read_namelist()

    def parse_system_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'SYSTEM' card.

        :return: a dictionary that stores the inputted information of 'SYSTEM' card
        """
        return SYSTEMNamelistParser(self.infile).read_namelist()

    def parse_electrons_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'ELECTRONS' card.

        :return: a dictionary that stores the inputted information of 'ELECTRONS' card
        """
        return ELECTRONSNamelistParser(self.infile).read_namelist()

    def parse_cell_parameters(self) -> np.ndarray:
        """
        Read 3 lines that follows 'CELL_PARAMETERS' string, so there must be no empty line between 'CELL_PARAMETERS' and
        the real cell parameters!

        :return: a numpy array that stores the cell parameters
        """
        cell_params = np.empty([3, 3])
        with open(self.infile, 'r') as f:
            for line in f:
                if 'CELL_PARAMETERS' in line.upper():
                    for i in range(3):
                        line = f.readline()
                        # if not line.strip():  # If there are blank lines
                        #     line = f.readline()
                        cell_params[i] = strs_to_floats(line.split())
        return cell_params

    def parse_k_points(self) -> Optional[KPoints]:
        """
        Find 'K_POINTS' line in the file, and read the k-mesh.
        If there is no 'K_POINTS' in file, it will read to end, and raise an error.

        We allow options and comments on the same line as 'K_POINTS':

        >>> test_strs = ['K_POINTS { crystal }','K_POINTS {crystal}','K_POINTS  crystal','K_POINTScrystal',\
        'K_POINTScrystal! This is a comment.','K_POINTS ! This is a comment.','K_POINTS']
        >>> [re.match("K_POINTS\s*{?\s*(\w*)\s*}?", s, re.IGNORECASE).groups()[0] for s in test_strs]
        ['crystal', 'crystal', 'crystal', 'crystal', 'crystal', '', '']

        :return: a named tuple defined above
        """
        option = 'tbipa'
        with open(self.infile, 'r') as f:
            for line in f:
                if 'K_POINTS' in line.upper():
                    option = re.match("K_POINTS\s*{?\s*(\w*)\s*}?", line, re.IGNORECASE).groups()[0]
                    sp: List[str] = f.readline().split()
                    grid = strs_to_ints(sp[0:3])
                    offsets = strs_to_ints(sp[3:7])
        return KPoints(grid=grid, offsets=offsets), option

    def parse_atomic_species(self):
        atmsp = []
        with open(self.infile, 'r') as f:
            for line in f:
                if 'ATOMIC_SPECIES' in line.upper():
                    for _ in range(int(self.parse_system_namelist()['ntyp'])):
                        line = f.readline()
                        if not line.strip():
                            line = f.readline()
                        name, mass, pseudopotential = line.strip().split()
                        atmsp.append(AtomicSpecies(name, float(mass), pseudopotential))
        return atmsp

    def parse_atomic_positions(self):
        atmpos = []
        option = 'alat'
        with open(self.infile, 'r') as f:
            for line in f:
                if 'ATOMIC_POSITIONS' in line.upper():
                    option = re.match("ATOMIC_POSITIONS\s*{?\s*(\w*)\s*}?", line, re.IGNORECASE).groups()[0]
                    for _ in range(int(self.parse_system_namelist()['nat'])):
                        atmpos.append(f.readline().strip().split())
        return atmpos, option

    def __str__(self):
        try:
            from beeprint import pp
            return str(pp(vars(self)))
        except ModuleNotFoundError:
            str(vars(self))

    __repr__ = __str__


# ====================================== The followings are output readers. ======================================
class PWscfOutputParser(SingleFileParser):
    def read_lattice_parameter(self) -> float:
        """
        Do not use it "as-is". This is a scaling number, not a 3x3 matrix containing the Cartesian coordinates.

        :return: a scaling number that multiplies `CELL_PARAMETERS` is the read Cartesian coordinates of the primitive
            vectors of a crystal.
        """
        return self._match_one_pattern("lattice parameter \(alat\)\s+=\s*(\d+\.\d+)", float)

    def read_total_energy(self) -> float:
        """
        Read total energy from the output file. The val unit is Rydberg.

        :return: total energy
        """
        return self._match_one_pattern("!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", float)

    def read_cell_volume(self) -> float:
        """
        Read the unit-cell volume, which is given at the beginning of the file.
        The val unit is $\text{bohr}^3$.

        :return: unit-cell volume
        """
        return self._match_one_pattern("unit-cell volume\s+=\s+(\d+\.\d+)", float)

    def read_pressure(self) -> float:
        """
        Read the pressure, which is at the bottom of the file.
        The val unit is kbar.

        :return: pressure at this volume
        """
        return self._match_one_pattern("P=\s+(-?\d+\.\d+)", float)

    def read_kinetic_energy_cutoff(self) -> float:
        """
        Read kinetic-energy cutoff, at the beginning of the file.

        :return: kinetic-energy cutoff, $E_\text{cut}$
        """
        return self._match_one_pattern("kinetic-energy cutoff\s+=\s*(-?\d+\.\d+)", float)

    def read_charge_density_cutoff(self) -> float:
        """
        Read charge density cutoff, at the beginning of the file.

        :return: charge density cutoff, $\rho_\text{cut}$
        """
        return self._match_one_pattern("charge density cutoff\s+=\s*(-?\d+\.\d+)", float)

    def read_atoms_num_per_cell(self) -> int:
        return self._match_one_pattern("number of atoms\/cell\s+=\s*(\d+)", int)

    def read_atoms_types_num(self) -> int:
        return self._match_one_pattern("number of atomic INPUTPH_types\s+=\s*(\d+)", int)

    def read_electrons_num(self) -> float:
        return self._match_one_pattern("number of electrons\s+=\s*(-?\d+\.\d+)", float)

    def read_mixing_beta(self) -> float:
        return self._match_one_pattern("mixing beta\s+=\s*(-?\d*\.\d+)", float)

    def read_nstep(self) -> int:
        return self._match_one_pattern("nstep\s+=\s*(\d+)", int)

    def read_iteration_num(self) -> int:
        """

        :return: the number of iterations used to reach self-consistent convergence
        """
        return self._match_one_pattern("number of iterations used\s+=\s*(\d+)", int)

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

        :param unit: There are 2 options, 'atomic' and 'kbar', where 'atomic' is the val one.
        :return: a 3x3 numpy array which contains the stress tensor
        """
        stress_atomic = np.zeros((3, 3))
        stress_kbar = np.zeros((3, 3))
        stress = {'atomic': stress_atomic, 'kbar': stress_kbar}
        with open(self.infile, 'r') as f:
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
            raise ValueError("Unknown coordinate system type! It can be either 'Cartesian' or 'crystal'!")

        with open(self.infile, 'r') as f, open(out_file, 'w') as g:
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

        with open(self.infile, 'r') as f:
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

    def read_processors_num(self) -> int:
        """
        Read how many processors were used in this calculation, which is given at the beginning of the file.

        :return: an integer denotes how many processors were used
        """
        return self._match_one_pattern("running on\s*(\d+)\s+", int)

    def read_symmetry_operations_num(self) -> int:
        return self._match_one_pattern("(\d+)\s+Sym\. Ops.*found", int)

    def read_fft_dimensions(self) -> List[int]:
        """

        :return: a list contains 3 integers that denotes the FFT grid we use
        """
        return list(map(int, self._match_one_pattern("FFT dimensions:\s+\(\s*(\d+),\s*(\d+),\s*(\d+)\)")))
