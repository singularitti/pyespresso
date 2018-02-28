#!/usr/bin/env python3

from functools import partialmethod

from pyespresso.core.namelist import CONTROLNamelistVariable, SYSTEMNamelistVariable, CELLNamelistVariable
from pyespresso.tools.strings import strings_to_integers


def partialclass(cls, *args, **kwds):
    class NewCls(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwds)

    return NewCls


PWSCF_GENERATOR_CONFIG = {
    'pseudopotential directory': partialclass(CONTROLNamelistVariable, 'pseudo_dir'),
    'prefix': partialclass(CONTROLNamelistVariable, 'prefix'),
    'fictitious cell mass': partialclass(CELLNamelistVariable, 'wmass'),
    'number of atoms': partialclass(SYSTEMNamelistVariable, 'nat'),
    'number of atomic types': partialclass(SYSTEMNamelistVariable, 'ntyp'),
    'energy cutoff': partialclass(SYSTEMNamelistVariable, 'ecutwfc'),
    'density cutoff': partialclass(SYSTEMNamelistVariable, 'ecutrho'),
    'smearing': partialclass(SYSTEMNamelistVariable, 'smearing'),
    'occupations': partialclass(SYSTEMNamelistVariable, 'occupations'),
    'degauss': partialclass(SYSTEMNamelistVariable, 'degauss'),
    'scratch folder': partialclass(CONTROLNamelistVariable, 'outdir'),
    'k-points': strings_to_integers,
    'shift': strings_to_integers
}


def pwscf_generator_config(raw, val):
    return PWSCF_GENERATOR_CONFIG[raw](val)
