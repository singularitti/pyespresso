#!/usr/bin/env python3

import operator
import os
import re
import warnings
from collections import namedtuple
from typing import *

import numpy as np
from lazy_property import *

from pyque.core.qe_input import AtomicSpecies, AtomicPosition, KPoints, PWscfInput
from pyque.lexer.simple import SimpleParser, NamelistParser
from pyque.meta.namelist import DEFAULT_CONTROL_NAMELIST, DEFAULT_SYSTEM_NAMELIST, DEFAULT_ELECTRONS_NAMELIST, \
    DEFAULT_IONS_NAMELIST, DEFAULT_CELL_NAMELIST
from pyque.meta.text import TextStream
from pyque.util.strings import strings_to_floats, strings_to_integers

# ========================================= What can be exported? =========================================
__all__ = ['SimpleParameter', 'to_text_file', 'PWscfInputParser', 'PWscfOutputParser']

# ================================= These are some type aliases or type definitions. =================================
SimpleParameter = NamedTuple('SimpleParameter', [('name', str), ('value', Union[str, int, float, bool]), ('type', str)])
RangeIndices = NamedTuple('RangeIndices', [('begin', int), ('end', int)])

# ========================================= define useful data structures =========================================
SimpleParameter: SimpleParameter = namedtuple('SimpleParameter', ['name', 'value', 'type'])
RangeIndices: RangeIndices = namedtuple('RangeIndices', ['begin', 'end'])


def to_text_file(obj: object, out_file: str):
    if isinstance(obj, PWscfInput):
        obj.to_text_file(out_file)
    else:
        raise TypeError('Input object is not an {0}!'.format('SCFStandardInput'))


# ====================================== The followings are input readers. ======================================
class CONTROLNamelistParser(NamelistParser):
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_CONTROL_NAMELIST)


class SYSTEMNamelistParser(NamelistParser):
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_SYSTEM_NAMELIST)


class ELECTRONSNamelistParser(NamelistParser):
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_ELECTRONS_NAMELIST)


class IONSNamelistParser(NamelistParser):
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_IONS_NAMELIST)


class CELLNamelistParser(NamelistParser):
    def __init__(self, infile):
        super().__init__(infile, DEFAULT_CELL_NAMELIST)


class PWscfInputParser(TextStream):
    """
    This class read an scf input file in, and parse it.
    """

    def __init__(self, instr: Optional[str] = None, infile: Optional[str] = None):
        self.linesep = "[\r\n]"
        self.namelist_sep = "/\s*[\r\n]"
        super().__init__(instr, infile)

    @LazyProperty
    def namelist_identifiers(self) -> List[str]:
        return ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL']

    @LazyProperty
    def card_identifiers(self) -> List[str]:
        return ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS', 'OCCUPATIONS', 'CONSTRAINTS',
                'ATOMIC_FORCES']

    @LazyProperty
    def namelist_identifier_positions(self, pos: int = 0, include_heading: bool = True,
                                      include_ending: bool = False) -> Dict[str, RangeIndices]:
        match_records = dict()
        identifiers = self.namelist_identifiers
        s = self.contents
        for pattern in identifiers:
            m0 = re.compile(pattern, flags=re.IGNORECASE).search(s, match_records.get(pattern, pos))
            if not m0:
                continue
            m1 = re.compile(self.namelist_sep, flags=re.IGNORECASE).search(s, m0.end())
            match_records[pattern] = RangeIndices(begin={True: m0.start(), False: m0.end()}[include_heading],
                                                  end={True: m1.end(), False: m1.start()}[include_ending])
        # The following one-line code first sort the ``dict`` *positions* by its values, i.e., a ``RangeIndices`` tuple.
        # Then we get a list of tuples, with first entries to be the identifiers and second indices to be
        # the indices.
        # For example, if ``x = {'a': 2, 'b': 3, 'c': 1}``, then
        # >>> sorted(x.items(), key=operator.itemgetter(1))
        # [('c', 1), ('a', 2), ('b', 3)]
        # Tuples are compared lexicographically using comparison of corresponding elements, thus compared with their
        # *begin* entry.
        m: List[Tuple[str, RangeIndices]] = sorted(match_records.items(), key=operator.itemgetter(1))
        return dict(m)

    @LazyProperty
    def card_identifier_positions(self):
        match_records = dict()
        identifiers = self.card_identifiers
        s = self.contents
        for pattern in identifiers:
            m = re.compile(pattern, flags=re.IGNORECASE).search(s, match_records.get(pattern, 500))
            if not m:
                continue
            match_records[pattern] = m.start()
        sorted_records = sorted(match_records.items(), key=operator.itemgetter(1))
        keys = list(map(operator.itemgetter(0), sorted_records))
        start_indices = list(map(operator.itemgetter(1), sorted_records))
        end_indices = list(map(lambda x: x - 1, start_indices[1:] + [len(s)]))
        return dict(zip(keys, (RangeIndices(begin=b, end=e) for b, e in zip(start_indices, end_indices))))

    def find_comment(self, include_heading: bool = True, include_ending: bool = False):
        pass
        # return self._section_with_bounds('!', '\R', include_heading, include_ending)

    def find_namelist(self, identifier):
        try:
            begin, end = self.namelist_identifier_positions[identifier]
        except KeyError:
            return ''
        return re.split(self.linesep, self.contents[begin:end])

    def find_card(self, identifier):
        try:
            begin, end = self.card_identifier_positions[identifier]
        except KeyError:
            return ''
        return re.split(self.linesep, self.contents[begin:end])

    def find_control_namelist(self):
        return self.find_namelist('&CONTROL')

    def find_system_namelist(self):
        return self.find_namelist('&SYSTEM')

    def find_electrons_namelist(self):
        return self.find_namelist('&ELECTRONS')

    def find_atomic_species(self):
        return self.find_card('ATOMIC_SPECIES')

    def find_atomic_positions(self):
        return self.find_card('ATOMIC_POSITIONS')

    def find_k_points(self):
        return self.find_card('K_POINTS')

    def find_cell_parameters(self):
        return self.find_card('CELL_PARAMETERS')

    def find_occupations(self):
        return self.find_card('OCCUPATIONS')

    def find_constraints(self):
        return self.find_card('CONSTRAINTS')

    def find_atomic_forces(self):
        return self.find_card('ATOMIC_FORCES')

    @LazyProperty
    def plain_control_namelist(self):
        return '\n'.join(self.find_control_namelist())

    @LazyProperty
    def plain_system_namelist(self):
        return '\n'.join(self.find_system_namelist())

    @LazyProperty
    def plain_electrons_namelist(self):
        return '\n'.join(self.find_electrons_namelist())

    @LazyProperty
    def plain_atomic_species(self):
        return '\n'.join(self.find_atomic_species())

    @LazyProperty
    def plain_atomic_positions(self):
        return '\n'.join(self.find_atomic_positions())

    @LazyProperty
    def plain_k_points(self):
        return '\n'.join(self.find_k_points())

    @LazyProperty
    def plain_cell_parameters(self):
        return '\n'.join(self.find_cell_parameters())

    @LazyProperty
    def plain_occupations(self):
        return '\n'.join(self.find_occupations())

    @LazyProperty
    def plain_constraints(self):
        return '\n'.join(self.find_constraints())

    @LazyProperty
    def plain_atomic_forces(self):
        return '\n'.join(self.find_atomic_forces())

    def parse_control_namelist(self) -> Dict[str, str]:
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

    def parse_atomic_species(self) -> Optional[List[AtomicSpecies]]:
        print(self.parse_system_namelist()['ntyp'])
        try:
            atom_types_number = self.parse_system_namelist()['ntyp'].value
        except KeyError:
            raise KeyError("The 'ntyp' parameter is not correctly given in SYSTEM namelist!")
        atomic_species = []
        generator: Iterator = self.stream_generator()
        for line in generator:
            if 'ATOMIC_SPECIES' in line.upper():
                for _ in range(atom_types_number):
                    if not line.strip():  # if this line is followed by an empty line
                        line = next(generator)
                    name, mass, pseudopotential = next(generator).strip().split()
                    atomic_species.append(AtomicSpecies(name, float(mass), pseudopotential))
                return atomic_species
        else:
            warnings.warn("No 'ATOMIC_SPECIES' is found in your input! Check it!")

    def parse_atomic_positions(self) -> Optional[Tuple[List[AtomicPosition], str]]:
        try:
            atoms_number = self.parse_system_namelist()['nat'].value
        except KeyError:
            raise KeyError("The 'nat' parameter is not correctly given in SYSTEM namelist!")
        atomic_positions = []
        generator: Iterator = self.stream_generator()
        for line in generator:
            if 'ATOMIC_POSITIONS' in line.upper():
                option = re.match("ATOMIC_POSITIONS\s*(?:\(|{)?\s*(\w*)\s*(?:\)|})?", line,
                                  re.IGNORECASE).group(1)
                if not option:  # if option is None:
                    warnings.warn("No option is found, default option 'alat' will be set!")
                    option = 'alat'
                for _ in range(atoms_number):
                    if not line.strip():  # if this line is followed by an empty line
                        line = next(generator)
                    name, coord1, coord2, coord3 = re.match(
                        "(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)",
                        next(generator).strip()).groups()
                    atomic_positions.append(
                        AtomicPosition(name, np.array(strings_to_floats([coord1, coord2, coord3]))))
                return atomic_positions, option
        else:
            warnings.warn("No 'ATOMIC_POSITIONS' is found in your input! Check it!")

    def parse_k_points(self) -> Optional[KPoints]:
        """
        Find 'K_POINTS' line in the file, and read the k-mesh.
        We allow options and comments on the same line as 'K_POINTS':

        >>> test_strs = ['K_POINTS { crystal }','K_POINTS {crystal}','K_POINTS  crystal','K_POINTScrystal', \
        'K_POINTScrystal! This is a comment.','K_POINTS ! This is a comment.','K_POINTS']
        >>> [re.match("K_POINTS\s*{?\s*(\w*)\s*}?", s, re.IGNORECASE).group(1) for s in test_strs]
        ['crystal', 'crystal', 'crystal', 'crystal', 'crystal', '', '']

        :return: a named tuple defined above
        """
        generator: Iterator = self.stream_generator()
        for line in generator:
            if 'K_POINTS' in line.upper():
                # The first parenthesized subgroup will be `option`.
                option = re.match("K_POINTS\s*(?:\(|{)?\s*(\w*)\s*(?:\)|})?", line, re.IGNORECASE).group(1)
                if not option:  # if option is None:
                    option = 'tpiba'
                try:
                    ks: List[int] = strings_to_integers(next(generator).split())
                except (ValueError, TypeError):
                    raise ValueError(
                        'This line is not a line of strings that can be converted into integers!')
                try:
                    grid, offsets = ks[0:3], ks[3:7]
                except IndexError:
                    raise IndexError('This line contains less than 6 integers!')
                return KPoints(grid=grid, offsets=offsets), option
        else:
            warnings.warn("No 'K_POINTS' is found in your input! Check it!")

    def parse_cell_parameters(self) -> Optional[Tuple[np.ndarray, str]]:
        """
        Read 3 lines that follows 'CELL_PARAMETERS' string, so there must be no empty line between 'CELL_PARAMETERS' and
        the real cell parameters!

        :return: a numpy array that stores the cell parameters
        """
        cell_params = np.empty([3, 3])
        generator: Iterator = self.stream_generator()
        for line in generator:
            if 'CELL_PARAMETERS' in line.upper():
                option = re.match("CELL_PARAMETERS\s*{?\s*(\w*)\s*}?", line, re.IGNORECASE).group(1)
                if not option:  # if option is None:
                    option = 'bohr'
                for i in range(3):
                    line = next(generator)
                    cell_params[i] = strings_to_floats(line.split())
                return cell_params, option
        else:
            warnings.warn(
                'Not specifying unit or lattice parameter is DEPRECATED and will no longer be allowed in the future!',
                category=DeprecationWarning)


# ====================================== The followings are output readers. ======================================
class PWscfOutputParser(SimpleParser):
    def read_lattice_parameter(self) -> float:
        """
        Do not use it "as-is". This is a scaling number, not a 3x3 matrix containing the Cartesian coordinates.

        :return: a scaling number that multiplies `CELL_PARAMETERS` is the read Cartesian coordinates of the primitive
            vectors of a crystal.
        """
        return self._match_one_pattern("lattice parameter \(alat\)\s+=\s*(\d+\.\d+)", float)

    def read_total_energy(self) -> float:
        """
        Read total energy from pyque.the output file. The val unit is Rydberg.

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
        Read the total stress, in 2 units, from pyque.the bottom of the output file.

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
        This method can be used to read how many k-points are involved in a PWscf calculation from pyque.a file
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
