#!/usr/bin/env python3
"""
:mod:`cards` -- 
========================================

.. module cards
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import re
from collections import namedtuple

import numpy as np

from pyespresso.tools.numbers import all_integer_like, all_float_like, is_float
from pyespresso.tools.strings import strings_to_integers, all_string_like, is_string_like, string_to_general_float, \
    strings_to_floats

# ========================================= What can be exported? =========================================
__all__ = ['AtomicSpecies', 'AtomicPosition', 'AutomaticKPoints', 'CellParameters']


class AtomicSpecies(namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])):
    """
    Note that the word 'species' serves as singular and plural both.
    So here though suffixed with a 's', it is for one atom and thus is singular."""

    def eval(self):
        mass = self.mass
        if is_string_like(self.mass):
            mass = string_to_general_float(self.mass)
        elif is_float(mass):
            pass
        else:
            raise TypeError
        return AtomicSpecies(self.name, mass, self.pseudopotential)

    def __str__(self) -> str:
        return '{0} {1} {2}'.format(self.name, self.mass, self.pseudopotential)

    __repr__ = __str__


class AutomaticKPoints(namedtuple('AutomaticKPoints', ['grid', 'offsets'])):
    __slots__ = ()

    def eval(self):
        grid, offsets = self.grid, self.offsets
        if all_string_like(grid):
            grid = strings_to_integers(grid)
        elif all_integer_like(grid):
            pass
        else:
            raise TypeError
        if all_string_like(offsets):
            offsets = strings_to_integers(offsets)
        elif all_integer_like(offsets):
            pass
        else:
            raise TypeError
        return AutomaticKPoints(grid, offsets)

    def __str__(self) -> str:
        return '{0} {1}'.format(' '.join(list(map(str, self.grid))), ' '.join(list(map(str, self.offsets))))

    __repr__ = __str__


class AtomicPosition(namedtuple('AtomicPosition', ['name', 'pos', 'if_pos'])):
    """
    *name*: label of the atom as specified in
    `ATOMIC_SPECIES <http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#ATOMIC_SPECIES>`_.

    *pos*: atomic positions, see `here <http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#idm1017>`_

    *if_pos*: component $i$ of the force for this atom is multiplied by *if_pos[i]*,
        which must be either $0$ or $1$.  Used to keep selected atoms and/or
        selected components fixed in MD dynamics or structural optimization run.
    """

    def eval(self):
        pos, if_pos = self.pos, self.if_pos
        if all_float_like(pos):
            pass
        elif all_string_like(pos):
            pos = strings_to_floats(pos)
        else:
            raise TypeError
        if all_integer_like(if_pos):
            pass
        elif all_string_like(if_pos):
            if_pos = strings_to_integers(if_pos)
        else:
            raise TypeError
        return AtomicPosition(self.name, pos, if_pos)

    def __str__(self) -> str:
        return '{0} {1} {2}'.format(self.name, ' '.join(list(map(str, self.pos))),
                                    ' '.join(list(map(str, self.if_pos))))

    __repr__ = __str__


class CellParameters:
    def __str__(self):
        return re.sub("[\[\]]", ' ', np.array2string(self,
                                                     formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))
