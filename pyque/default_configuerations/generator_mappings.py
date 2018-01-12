#!/usr/bin/env python3

from functools import partialmethod

from pyque.meta.parameter import CONTROLNamelistParameter, SYSTEMNamelistParameter, CELLNamelistParameter
from pyque.miscellaneous.strings import strs_to_ints


def partialclass(cls, *args, **kwds):
    class NewCls(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwds)

    return NewCls


PWSCF_GENERATOR_CONFIG = {
    'pseudopotential directory': partialclass(CONTROLNamelistParameter, 'pseudo_dir'),
    'prefix': partialclass(CONTROLNamelistParameter, 'prefix'),
    'fictitious cell mass': partialclass(CELLNamelistParameter, 'wmass'),
    'number of atoms': partialclass(SYSTEMNamelistParameter, 'nat'),
    'number of atomic types': partialclass(SYSTEMNamelistParameter, 'ntyp'),
    'energy cutoff': partialclass(SYSTEMNamelistParameter, 'ecutwfc'),
    'density cutoff': partialclass(SYSTEMNamelistParameter, 'ecutrho'),
    'smearing': partialclass(SYSTEMNamelistParameter, 'smearing'),
    'occupations': partialclass(SYSTEMNamelistParameter, 'occupations'),
    'degauss': partialclass(SYSTEMNamelistParameter, 'degauss'),
    'scratch folder': partialclass(CONTROLNamelistParameter, 'outdir'),
    'k-points': strs_to_ints,
    'shift': strs_to_ints
}


def pwscf_generator_config(raw, val):
    return PWSCF_GENERATOR_CONFIG[raw](val)
