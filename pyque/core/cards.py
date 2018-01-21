#!/usr/bin/env python3
"""
:mod:`cards` -- 
========================================

.. module cards
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from collections import namedtuple

# ========================================= What can be exported? =========================================
__all__ = ['AtomicSpecies', 'AtomicPosition', 'KPoints']


class AtomicSpecies(namedtuple('AtomicSpecies', ['name', 'mass', 'pseudopotential'])):
    __doc__ = """Note that the word 'species' serves as singular and plural both.
    So here though suffixed with a 's', it is for one atom and thus is singular."""


class KPoints(namedtuple('KPoints', ['grid', 'offsets'])):
    pass


class AtomicPosition(namedtuple('AtomicPosition', ['name', 'x', 'y', 'z', 'if_pos1', 'if_pos2', 'if_pos3'])):
    pass
