#!/usr/bin/env python3
# created on Jan 6, 2018 at 02:33 by Qi Zhang

from functools import partialmethod

from data_models.parameters import *
from miscellaneous.string import strs_to_ints


def partialclass(cls, *args, **kwds):
    class NewCls(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwds)

    return NewCls


PWSCF_GENERATOR_CONFIG = {
    'pseudopotential directory': partialclass(CONTROLParameter, 'pseudo_dir'),
    'prefix': partialclass(CONTROLParameter, 'prefix'),
    'fictitious cell mass': partialclass(CELLParameter, 'wmass'),
    'number of atoms': partialclass(SYSTEMParameter, 'nat'),
    'number of atomic types': partialclass(SYSTEMParameter, 'ntyp'),
    'energy cutoff': partialclass(SYSTEMParameter, 'ecutwfc'),
    'density cutoff': partialclass(SYSTEMParameter, 'ecutrho'),
    'smearing': partialclass(SYSTEMParameter, 'smearing'),
    'occupations': partialclass(SYSTEMParameter, 'occupations'),
    'degauss': partialclass(SYSTEMParameter, 'degauss'),
    'scratch folder': partialclass(CONTROLParameter, 'outdir'),
    'k-points': strs_to_ints,
    'shift': strs_to_ints
}


def pwscf_generator_config(raw, val):
    return PWSCF_GENERATOR_CONFIG[raw](val)
