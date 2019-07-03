#!/usr/bin/env python3
# Created at Mar 5, 2018, by Qi Zhang

import operator
import re
import warnings
from collections import namedtuple, OrderedDict
from typing import *

import numpy as np
from f90nml.tokenizer import Tokenizer as NamelistLexer
from lazy_property import LazyProperty

from pyespresso.core.cards import AtomicSpecies, AtomicPosition, AutomaticKPoints
from pyespresso.meta.text import TextStream


# ========================================= define useful data structures =========================================
class RangeIndices(namedtuple('RangeIndices', ['begin', 'end'])):
    def __str__(self) -> str:
        return "'begin: {0}, end: {1}'".format(self.begin, self.end)


# ====================================== The followings are core data structures. ======================================
class PWscfInputLexer:
    """
    This class reads a standard Quantum ESPRESSO PWscf input file or string in, and lex it.
    """

    def __init__(self, inp: Optional[str] = None, **kwargs):
        self.newline = "[\r\n,]"  # TODO: This will fail when ',' is inside a value of a parameter.
        self.namelist_sep = "/\s*[\r\n]"
        self.__text_stream = TextStream(inp, **kwargs)

    @property
    def namelist_identifiers(self) -> List[str]:
        return ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL']

    @property
    def card_identifiers(self) -> List[str]:
        return ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS', 'OCCUPATIONS', 'CONSTRAINTS',
                'ATOMIC_FORCES']

    def __get_card_identifier_positions(self, pos: Optional[int] = None) -> MutableMapping[str, RangeIndices]:
        """


        :param pos:
        :return:
        """
        match_records = dict()
        identifiers: List[str] = self.card_identifiers
        s: str = self.__text_stream.content
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

    def __get_card(self, identifier) -> Optional[List[str]]:
        if identifier in self.cards_found:
            begin, end = self.__get_card_identifier_positions()[identifier]
            return re.split(self.newline, self.__text_stream.content[begin:end])
        else:
            warnings.warn("Identifier '{0}' is not found in input!".format(identifier), stacklevel=2)

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

    def lex_namelist(self):
        """
        A generic method to read a namelist.
        Note you cannot write more than one parameter in each line!

        :return: a dictionary that stores the inputted information of the intended card
        """
        return NamelistLexer(self.__text_stream)

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
