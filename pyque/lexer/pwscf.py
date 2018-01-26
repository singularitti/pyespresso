#!/usr/bin/env python3

import operator
import os
import re
import warnings
from collections import namedtuple, OrderedDict
from typing import *

import numpy as np
from lazy_property import LazyProperty

from pyque.core.cards import AtomicSpecies, AtomicPosition, AutomaticKPoints, CellParameters
from pyque.meta.text import TextStream

# ========================================= What can be exported? =========================================
__all__ = ['PWscfInputLexer', 'PWscfOutputLexer']


# ========================================= define useful data structures =========================================
class RangeIndices(namedtuple('RangeIndices', ['begin', 'end'])):
    def __str__(self) -> str:
        return "'begin: {0}, end: {1}'".format(self.begin, self.end)


# ========================================= define useful functions =========================================


# ====================================== The followings are core data structures. ======================================
class PWscfInputLexer:
    """
    This class reads a standard Quantum ESPRESSO PWscf input file or string in, and lex it.
    """

    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None):
        self.linesep = "[\r\n,]"  # TODO: This will fail when ',' is inside a value of a parameter.
        self.namelist_sep = "/\s*[\r\n]"
        self.__text_stream = TextStream(inp=instream, infile=infile)

    @property
    def namelist_identifiers(self) -> List[str]:
        return ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL']

    @property
    def card_identifiers(self) -> List[str]:
        return ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS', 'OCCUPATIONS', 'CONSTRAINTS',
                'ATOMIC_FORCES']

    def __get_namelist_identifier_positions(self, pos: int = 0, include_heading: bool = True,
                                            include_ending: bool = False) -> MutableMapping[str, RangeIndices]:
        """
        For example, a typical returned result will look like

        .. code-block:: python

            {'&CONTROL': 'begin: 1, end: 440', '&SYSTEM': 'begin: 443, end: 629', '&ELECTRONS': 'begin: 632, end: 684'}

        :param pos:
        :param include_heading:
        :param include_ending:
        :return:
        """
        match_records = dict()
        identifiers: List[str] = self.namelist_identifiers
        s: str = self.__text_stream.contents
        for pattern in identifiers:
            # ``re.compile`` will produce a regular expression object, on which we can use its ``search`` method.
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
        return OrderedDict(m)

    def __get_card_identifier_positions(self, pos: Optional[int] = None) -> MutableMapping[str, RangeIndices]:
        """


        :param pos:
        :return:
        """
        match_records = dict()
        identifiers: List[str] = self.card_identifiers
        s: str = self.__text_stream.contents
        if not pos:
            pos = list(self.__get_namelist_identifier_positions().values())[-1].end
        for pattern in identifiers:
            m = re.compile(pattern, flags=re.IGNORECASE).search(s, match_records.get(pattern, pos))
            if not m:
                continue
            match_records[pattern] = m.start()
        sorted_records = sorted(match_records.items(), key=operator.itemgetter(1))
        keys = list(map(operator.itemgetter(0), sorted_records))
        start_indices = list(map(operator.itemgetter(1), sorted_records))
        end_indices = list(map(lambda x: x - 1, start_indices[1:] + [len(s)]))
        return OrderedDict(zip(keys, (RangeIndices(begin=b, end=e) for b, e in zip(start_indices, end_indices))))

    @LazyProperty
    def namelists_found(self) -> Optional[Set[str]]:
        """
        Check whether an input contains all the namelists necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all namelists found in the input will be returned.

        :return: All namelists found in the input.
        """
        keys = list(self.__get_namelist_identifier_positions().keys())
        for i, k in enumerate(keys):
            if not k == self.namelist_identifiers[i]:
                # Namelists must appear in the order given below.
                raise RuntimeError("Namelists Must be in order 'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'!")
        else:
            return set(map(lambda s: s[1:], keys))

    @LazyProperty
    def cards_found(self) -> Optional[Set[str]]:
        """
        Check whether an input contains all the cards necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all cards found in the input will be returned.

        :return: All cards found in the input.
        """
        keys = set(self.__get_card_identifier_positions().keys())
        if keys < {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'}:
            # 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', and 'K_POINTS' are necessary cards.
            warnings.warn('Not enough necessary namelists given!')
        else:
            return keys

    def get_comment(self, include_heading: bool = True, include_ending: bool = False):
        pass
        # return self._section_with_bounds('!', '\R', include_heading, include_ending)

    def __get_namelist(self, group_name) -> Optional[List[str]]:
        # TODO: This will fail when ',' is inside a value of a parameter.
        identifier = '&' + group_name
        if group_name in self.namelists_found:
            begin, end = self.__get_namelist_identifier_positions()[identifier]
            return re.split(self.linesep, self.__text_stream.contents[begin:end])
        else:
            warnings.warn("Identifier '{0}' is not found in input!".format(identifier), stacklevel=2)

    def __get_card(self, identifier) -> Optional[List[str]]:
        if identifier in self.cards_found:
            begin, end = self.__get_card_identifier_positions()[identifier]
            return re.split(self.linesep, self.__text_stream.contents[begin:end])
        else:
            warnings.warn("Identifier '{0}' is not found in input!".format(identifier), stacklevel=2)

    def get_control_namelist(self):
        return self.__get_namelist('CONTROL')

    def get_system_namelist(self):
        return self.__get_namelist('SYSTEM')

    def get_electrons_namelist(self):
        return self.__get_namelist('ELECTRONS')

    def get_ions_namelist(self):
        return self.__get_namelist('IONS')

    def get_cell_namelist(self):
        return self.__get_namelist('CELL')

    def get_atomic_species(self):
        return self.__get_card('ATOMIC_SPECIES')

    def get_atomic_positions(self):
        return self.__get_card('ATOMIC_POSITIONS')

    def get_k_points(self):
        return self.__get_card('K_POINTS')

    def get_cell_parameters(self):
        return self.__get_card('CELL_PARAMETERS')

    def get_occupations(self):
        return self.__get_card('OCCUPATIONS')

    def get_constraints(self):
        return self.__get_card('CONSTRAINTS')

    def get_atomic_forces(self):
        return self.__get_card('ATOMIC_FORCES')

    def lex_namelist(self, group_name: str) -> Optional[Dict[str, str]]:
        """
        A generic method to read a namelist.
        Note you cannot write more than one parameter in each line!

        :return: a dictionary that stores the inputted information of the intended card
        """
        result = dict()
        for line in self.__get_namelist(group_name):  # Read each line in the namelist until '/'
            s: str = line.strip()
            # Use '=' as the delimiter, split the stripped line into a key and a value.
            # Skip this line if a line starts with '&' (namelist caption) or '!' (comment) or
            # this line is empty ('').
            if s.startswith('&') or s.startswith('!') or not s:
                continue
            k, v = s.split('=', maxsplit=1)
            k: str = k.strip()
            v: str = v.strip().rstrip(',').strip().split('!')[0]  # Ignore trailing comma of the line
            result.update({k: v.strip()})
        return result

    def lex_control_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'CONTROL' namelist.

        :return: A dictionary that stores the information of 'CONTROL' namelist.
        """
        return self.lex_namelist('&CONTROL')

    def lex_system_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'SYSTEM' namelist.

        :return: A dictionary that stores the inputted information of 'SYSTEM' namelist.
        """
        return self.lex_namelist('&SYSTEM')

    def lex_electrons_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'ELECTRONS' namelist.

        :return: A dictionary that stores the information of 'ELECTRONS' namelist.
        """
        return self.lex_namelist('&ELECTRONS')

    def lex_ions_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'IONS' namelist.

        :return: A dictionary that stores the information of 'IONS' namelist.
        """
        return self.lex_namelist('&IONS')

    def lex_cell_namelist(self) -> Dict[str, str]:
        """
        Read everything that falls within 'CELL' namelist.

        :return: A dictionary that stores the information of 'IONS' namelist.
        """
        return self.lex_namelist('&CELL')

    def lex_atomic_species(self) -> Optional[List[AtomicSpecies]]:
        s: Optional[List[str]] = self.get_atomic_species()
        if not s:  # If the returned result is ``None``.
            warnings.warn("'ATOMIC_SPECIES' not found in input!", stacklevel=2)
        else:
            atomic_species = []
            for line in s[1:]:
                # Skip the title line, any empty line, or a line of comment.
                if not line.strip() or line.strip().startswith('!'):
                    continue
                match = re.match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", line.strip())
                if match is None:
                    warnings.warn("No match found in the line {0}!".format(line))
                else:
                    name, mass, pseudopotential = match.groups()
                    atomic_species.append(AtomicSpecies(name, mass, pseudopotential))
            return atomic_species

    def lex_atomic_positions(self) -> Optional[Tuple[List[AtomicPosition], str]]:
        s: Optional[List[str]] = self.get_atomic_positions()
        if not s:  # If the returned result is ``None``.
            warnings.warn("'ATOMIC_POSITIONS' not found in input!", stacklevel=2)
        else:
            atomic_positions = []
            title_line = s[0]
            match = re.match("ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?", title_line, flags=re.IGNORECASE)
            if match is None:
                raise RuntimeError("No match found in the line '{0}'! Something went wrong!".format(title_line))
            option = match.group(1)
            if option == '':
                warnings.warn("No option is found, default option 'alat' will be set! "
                              "Not specifying units is DEPRECATED and will no longer be allowed in the future",
                              category=DeprecationWarning)
                option = 'alat'
            for line in s[1:]:
                # If this line is an empty line or a line of comment.
                if line.strip() == '' or line.strip().startswith('!'):
                    continue
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                if re.match("{.*}", line):
                    match = re.match(
                        "(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*{\s*([01])?\s*([01])?\s*([01])?\s*}",
                        line.strip())
                    name, x, y, z, if_pos1, if_pos2, if_pos3 = match.groups()
                    atomic_positions.append(AtomicPosition(name, [x, y, z], [if_pos1, if_pos2, if_pos3]))
                else:
                    match = re.match("(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", line.strip())
                    if match is None:
                        warnings.warn("No match found in the line {0}!".format(line))
                    else:
                        name, x, y, z = match.groups()
                        atomic_positions.append(AtomicPosition(name, [x, y, z], ['1', '1', '1']))
            return atomic_positions, option

    # TODO: finish this method
    def lex_k_points(self) -> Union[None, str, Tuple[AutomaticKPoints, str]]:
        """
        Find 'K_POINTS' line in the file, and read the k-mesh.
        We allow options and comments on the same line as 'K_POINTS':

        >>> test_strs = ['K_POINTS { crystal }','K_POINTS {crystal}','K_POINTS  crystal','K_POINTScrystal', \
        'K_POINTScrystal! This is a comment.','K_POINTS ! This is a comment.','K_POINTS']
        >>> [re.match("K_POINTS\s*{?\s*(\w*)\s*}?", s, re.IGNORECASE).group(1) for s in test_strs]
        ['crystal', 'crystal', 'crystal', 'crystal', 'crystal', '', '']

        :return: a named tuple defined above
        """
        s: Optional[List[str]] = self.get_k_points()
        title_line = s[0]
        match = re.match("K_POINTS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?", title_line, flags=re.IGNORECASE)
        if match is None:
            raise RuntimeError("Match not found! Check your option!")
        option = match.group(1)  # The first parenthesized subgroup will be `option`.
        if option == '':
            raise RuntimeError("Option is not given! you must give one!")
        elif option == 'gamma':
            return option
        elif option == 'automatic':
            for line in s[1:]:
                if line.strip() == '' or line.strip().startswith('!'):
                    continue
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                line = line.split()
                grid, offsets = line[0:3], line[3:7]
                return AutomaticKPoints(grid=grid, offsets=offsets), option
        elif option in {'tpiba', 'crystal', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c'}:
            return NotImplemented
        else:
            raise ValueError("Unknown option '{0}' given!".format(option))

    def lex_cell_parameters(self):
        """
        Read 3 lines that follows 'CELL_PARAMETERS' string, so there must be no empty line between 'CELL_PARAMETERS' and
        the real cell parameters!

        :return: a numpy array that stores the cell parameters
        """
        if not self.get_cell_parameters():  # If returned result is ``None``.
            warnings.warn("'CELL_PARAMETERS' not found in input!", stacklevel=2)
        else:
            cell_params = []
            title_line = self.get_cell_parameters()[0]
            match = re.match("CELL_PARAMETERS\s*{?\s*(\w*)\s*}?", title_line, re.IGNORECASE)
            if match is None:
                # The first line should be 'CELL_PARAMETERS blahblahblah', if it is not, either the regular expression
                # wrong or something worse happened.
                raise RuntimeError("No match found! Check you 'CELL_PARAMETERS' line!")
            option = match.group(1)
            if option == '':
                warnings.warn('Not specifying unit is DEPRECATED and will no longer be allowed in the future!',
                              category=DeprecationWarning)
                option = 'bohr'
            for line in self.get_cell_parameters()[1:]:
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                if re.match("(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", line.strip()):
                    v1, v2, v3 = re.match("(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", line.strip()).groups()
                    cell_params.append([v1, v2, v3])
            return np.array(cell_params), option


# ====================================== The followings are output readers. ======================================
class PWscfOutputLexer:
    def __init__(self, instream: Optional[str] = None, infile: Optional[str] = None, encoding: Optional[str] = None,
                 newline='\n'):
        self.text_stream = TextStream(instream, infile, encoding, newline)

    def _match_one_pattern(self, pattern: str, *args, **kwargs) -> Union[None, Any, List]:
        """
        This method matches a pattern which exists only once in the file.

        :param pattern: a regular expression that you want to match
        :param args: a wrapper function which determines the returned type of value
        :return: Determined by the `wrapper`, the value you want to grep out from pyque.the file.
        """
        s = self.text_stream.contents()
        match: Optional[List[str]] = re.findall(pattern, s,
                                                **kwargs)  # `match` is either an empty list or a list of strings.
        if match:
            if len(args) == 0:  # If no wrapper argument is given, return directly the matched string
                return match
            elif len(args) == 1:  # If wrapper argument is given, i.e., not empty, then apply wrapper to the match
                wrapper, = args
                return [wrapper(m) for m in match]
            else:
                raise TypeError('Multiple wrappers are given! Only one should be given!')
        else:  # If no match is found
            print('Pattern {0} not found in string {1}!'.format(pattern, s))

    def lex_lattice_parameter(self) -> float:
        """
        Do not use it "as-is". This is a scaling number, not a 3x3 matrix containing the Cartesian coordinates.

        :return: a scaling number that multiplies `CELL_PARAMETERS` is the read Cartesian coordinates of the primitive
            vectors of a crystal.
        """
        return self._match_one_pattern("lattice parameter \(alat\)\s+=\s*(\d+\.\d+)", float)

    def lex_total_energy(self) -> float:
        """
        Read total energy from pyque.the output file. The val unit is Rydberg.

        :return: total energy
        """
        return self._match_one_pattern("!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", float)

    def lex_cell_volume(self) -> float:
        """
        Read the unit-cell volume, which is given at the beginning of the file.
        The val unit is $\text{bohr}^3$.

        :return: unit-cell volume
        """
        return self._match_one_pattern("unit-cell volume\s+=\s+(\d+\.\d+)", float)

    def lex_pressure(self) -> float:
        """
        Read the pressure, which is at the bottom of the file.
        The val unit is kbar.

        :return: pressure at this volume
        """
        return self._match_one_pattern("P=\s+(-?\d+\.\d+)", float)

    def lex_kinetic_energy_cutoff(self) -> float:
        """
        Read kinetic-energy cutoff, at the beginning of the file.

        :return: kinetic-energy cutoff, $E_\text{cut}$
        """
        return self._match_one_pattern("kinetic-energy cutoff\s+=\s*(-?\d+\.\d+)", float)

    def lex_charge_density_cutoff(self) -> float:
        """
        Read charge density cutoff, at the beginning of the file.

        :return: charge density cutoff, $\rho_\text{cut}$
        """
        return self._match_one_pattern("charge density cutoff\s+=\s*(-?\d+\.\d+)", float)

    def lex_atoms_num_per_cell(self) -> int:
        return self._match_one_pattern("number of atoms\/cell\s+=\s*(\d+)", int)

    def lex_atoms_types_num(self) -> int:
        return self._match_one_pattern("number of atomic INPUTPH_types\s+=\s*(\d+)", int)

    def lex_electrons_num(self) -> float:
        return self._match_one_pattern("number of electrons\s+=\s*(-?\d+\.\d+)", float)

    def lex_mixing_beta(self) -> float:
        return self._match_one_pattern("mixing beta\s+=\s*(-?\d*\.\d+)", float)

    def lex_nstep(self) -> int:
        return self._match_one_pattern("nstep\s+=\s*(\d+)", int)

    def lex_iteration_num(self) -> int:
        """

        :return: the number of iterations used to reach self-consistent convergence
        """
        return self._match_one_pattern("number of iterations used\s+=\s*(\d+)", int)

    def lex_kinetic_stress(self) -> np.ndarray:
        return self._lex_stress('kinetic')

    def lex_local_stress(self) -> np.ndarray:
        return self._lex_stress('local')

    def lex_nonlocal_stress(self) -> np.ndarray:
        return self._lex_stress('nonloc.')

    def lex_hartree_stress(self) -> np.ndarray:
        return self._lex_stress('hartree')

    def lex_exc_cor_stress(self) -> np.ndarray:
        return self._lex_stress('exc-cor')

    def lex_corecor_stress(self) -> np.ndarray:
        return self._lex_stress('corecor')

    def lex_ewald_stress(self) -> np.ndarray:
        return self._lex_stress('ewald')

    def lex_hubbard_stress(self) -> np.ndarray:
        return self._lex_stress('hubbard')

    def lex_london_stress(self) -> np.ndarray:
        return self._lex_stress('london')

    def lex_xdm_stress(self) -> np.ndarray:
        return self._lex_stress('XDM')

    def lex_dft_nl_stress(self) -> np.ndarray:
        return self._lex_stress('dft-nl')

    def lex_ts_vdw_stress(self) -> np.ndarray:
        return self._lex_stress('TS-vdW')

    def lex_total_stress(self, unit: Optional[str] = 'atomic') -> np.ndarray:
        """
        Read the total stress, in 2 units, from pyque.the bottom of the output file.

        :param unit: There are 2 options, 'atomic' and 'kbar', where 'atomic' is the val one.
        :return: a 3x3 numpy array which contains the stress tensor
        """
        stress_atomic = np.zeros((3, 3))
        stress_kbar = np.zeros((3, 3))
        stress = {'atomic': stress_atomic, 'kbar': stress_kbar}
        for line in self.text_stream.stream_generator():
            if 'total   stress' in line:
                for i in range(3):  # Read a 3x3 matrix
                    line = next(line)
                    stress_atomic[i][:] = list(map(float, line.split()))[0:3]
                    stress_kbar[i][:] = list(map(float, line.split()))[3:6]
        return stress[unit]

    def lex_k_coordinates(self, out_file: str, coordinate_system: Optional[str] = 'crystal'):
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

        with open(out_file, 'w') as g:
            flag = False  # If flag is true, read line and match pattern, if not true, no need to match pattern
            k_count = 0  # Count how many k-points have been read
            k_num = None  # How many k-points in total, given by Quantum ESPRESSO
            for line in self.text_stream.stream_generator():
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

    def _lex_stress(self, name) -> np.ndarray:
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

        for line in self.text_stream.stream_generator():
            if name in line and 'stress' in line:  # Read the first line of the matrix
                try:
                    stress[0][:] = list(map(float, re.findall(name + reg1, line)[0][:]))
                except IndexError:
                    raise IndexError(error_message)
                try:
                    for i in range(1, 3):  # Read the rest 2 lines of matrix
                        line = next(line)
                        stress[i][:] = list(map(float, re.findall(reg2, line)[0][:]))
                except IndexError:
                    raise IndexError(error_message)
        return stress

    def lex_processors_num(self) -> int:
        """
        Read how many processors were used in this calculation, which is given at the beginning of the file.

        :return: an integer denotes how many processors were used
        """
        return self._match_one_pattern("running on\s*(\d+)\s+", int)

    def lex_symmetry_operations_num(self) -> int:
        return self._match_one_pattern("(\d+)\s+Sym\. Ops.*found", int)

    def lex_fft_dimensions(self) -> List[int]:
        """

        :return: a list contains 3 integers that denotes the FFT grid we use
        """
        return list(map(int, self._match_one_pattern("FFT dimensions:\s+\(\s*(\d+),\s*(\d+),\s*(\d+)\)")))
