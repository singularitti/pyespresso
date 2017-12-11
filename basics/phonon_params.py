#!/usr/bin/env python3
# created at Nov 20, 2017 3:17 PM by Qi Zhang

from collections import defaultdict

from basics.parameter import *

# Type alias
_INPUTPHDefaultParameter = DefaultDict[str, Union[str, float, int, bool]]

names = ['amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph', 'alpha_mix', 'nmix_ph', 'verbosity', 'reduce_io',
         'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil', 'lrpa', 'lnoloc', 'trans', 'lraman', 'eth_rps',
         'eth_ns', 'dek', 'recover', 'low_directory_check', 'only_init', 'qplot', 'q2d', 'q_in_band_form',
         'electron_phonon', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol', 'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir',
         'search_sym', 'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3', 'start_irr', 'last_irr',
         'nat_todo', 'modenum', 'start_q', 'last_q', 'dvscf_star', 'drho_star']

default_values = [0, './', 'pwscf', 100, 1e-12, 0.7, 4, 'val', False,
                  1e7, 'matdyn', '', '', False, False, False, True, False, 1e-9,
                  1e-12, 1e-3, False, False, False, False, False, False,
                  '', False, False, False, False, False, False, False, False, False, False,
                  True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3,
                  0, 0, 1, 3, 'disabled', 'disabled']

_INPUTPHDefaultParameter: _INPUTPHDefaultParameter = defaultdict(tuple)
for name, default_value in zip(names, default_values):
    _INPUTPHDefaultParameter.update({name: default_value})


class INPUTPHParameter(Parameter):
    def __init__(self, name, value):
        super().__init__(name, value)
        self._default_value = _INPUTPHDefaultParameter[name]

    @property
    def default_value(self):
        """
        The default value given by Quantum ESPRESSO. Read-only property.

        :return: the default value given by Quantum ESPRESSO
        """
        return self._default_value

# if __name__ == '__main__':
#     eg = PhononParam('amass', '56')
#     print(type(eg.value) == str)
#     print(type(eg.value) == float)
