#!/usr/bin/env python3
# created at Nov 20, 2017 3:17 PM by Qi Zhang

from collections import namedtuple

from miscellaneous.param import *

# It has to be a list but not a set, or there will be alignment problem.
INPUTPH_card = [
    'amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph', 'alpha_mix(niter)', 'nmix_ph', 'verbosity', 'reduce_io',
    'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil', 'lrpa', 'lnoloc', 'trans', 'lraman', 'eth_rps',
    'eth_ns', 'dek', 'recover', 'low_directory_check', 'only_init', 'qplot', 'q2d', 'q_in_band_form',
    'electron_phonon', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol', 'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir',
    'search_sym', 'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3', 'start_irr', 'last_irr',
    'nat_todo', 'modenum', 'start_q', 'last_q', 'dvscf_star', 'drho_star'
]

# It has to be a list but not a set, or same strings will be removed.
# Default card types for corresponding name
INPUTPH_types = [
    'float', 'str', 'str', 'int', 'float', 'float', 'int', 'str', 'bool',
    'float', 'str', 'str', 'str', 'bool', 'bool', 'bool', 'bool', 'bool', 'float',
    'float', 'float', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool',
    'str', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool', 'bool',
    'bool', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int',
    'int', 'int', 'int', 'int', 'str', 'str'
]

# Default values for corresponding name and type
INPUTPH_default_values = [
    0, './', 'pwscf', 100, 1e-12, 0.7, 4, 'val', False,
    1e7, 'matdyn', '', '', False, False, False, True, False, 1e-9,
    1e-12, 1e-3, False, False, False, False, False, False,
    '', False, False, False, False, False, False, False, False, False, False,
    True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3,
    0, 0, 1, 3, 'disabled', 'disabled'
]

# Define a namedtuple to store all default values indicated by Quantum ESPRESSO
INPUTPH_item: NamedTuple = namedtuple('INPUTPH_item', ['name', 'default_value', 'type'])

# Form a dictionary of default values
INPUTPH_namedtuples: Dict[str, NamedTuple] = {name: INPUTPH_item(name, val, typ) for name, val, typ in
                                              zip(INPUTPH_card, INPUTPH_default_values, INPUTPH_types)}


class PhononParam(Param):
    def __init__(self, name: str, raw_val: str):
        """
        Generate an INPUTPH object, which stores user given name, and value.

        :param name: the name given by user in the INPUTPH card
        :param raw_val: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        super().__init__(name, raw_val, INPUTPH_namedtuples)


if __name__ == '__main__':
    eg = PhononParam('amass', '56')
    print(type(eg.raw_value) == str)
    print(type(eg.value) == float)
