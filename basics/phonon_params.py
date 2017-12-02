#!/usr/bin/env python3
# created at Nov 20, 2017 3:17 PM by Qi Zhang

from collections import namedtuple

from basics.param import *

# Define a namedtuple to store all default values indicated by Quantum ESPRESSO
INPUTPH_param: NamedTuple = namedtuple('INPUTPH_param', ['name', 'default_value', 'type'])


class INPUTPHNamelist:
    @property
    def names(self):
        # It has to be a list but not a set, or there will be alignment problem.
        return ['amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph', 'alpha_mix', 'nmix_ph', 'verbosity', 'reduce_io',
                'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil', 'lrpa', 'lnoloc', 'trans', 'lraman', 'eth_rps',
                'eth_ns', 'dek', 'recover', 'low_directory_check', 'only_init', 'qplot', 'q2d', 'q_in_band_form',
                'electron_phonon', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol', 'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir',
                'search_sym', 'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3', 'start_irr', 'last_irr',
                'nat_todo', 'modenum', 'start_q', 'last_q', 'dvscf_star', 'drho_star']

    @property
    def param_types(self) -> List[type]:
        """
        Default namelist's types for corresponding name

        :return:
        """
        # It has to be a list but not a set, or same strings will be removed.
        return [float, str, str, int, float, float, int, str, bool,
                float, str, str, str, bool, bool, bool, bool, bool, float,
                float, float, bool, bool, bool, bool, bool, bool,
                str, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool,
                bool, int, int, int, int, int, int, int, int, int, int, int,
                int, int, int, int, str, str]

    @property
    def default_values(self):
        """
        Default values for corresponding name and type

        :return:
        """
        return [0, './', 'pwscf', 100, 1e-12, 0.7, 4, 'val', False,
                1e7, 'matdyn', '', '', False, False, False, True, False, 1e-9,
                1e-12, 1e-3, False, False, False, False, False, False,
                '', False, False, False, False, False, False, False, False, False, False,
                True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3,
                0, 0, 1, 3, 'disabled', 'disabled']


class INPUTPHDefaultParams:
    @property
    def params(self) -> Dict[str, NamedTuple]:
        d = {}
        ipn = INPUTPHNamelist()
        names, default_values, param_types = ipn.names, ipn.default_values, ipn.param_types
        for name, val, typ in zip(names, default_values, param_types):
            d.update({name: INPUTPH_param(name, val, typ)})
        return d


print(INPUTPHDefaultParams().params)


class INPUTPHParam(Param):
    @property
    def default_value(self):
        """
        The default value given by Quantum ESPRESSO. Read-only property.

        :return: the default value given by Quantum ESPRESSO
        """
        return self._default_value

    @property
    def type(self) -> type:
        """
        The default type defined by Quantum ESPRESSO. Read-only property.

        :return: a string given by Quantum ESPRESSO manual
        """
        return self._type

# if __name__ == '__main__':
#     eg = PhononParam('amass', '56')
#     print(type(eg.raw_value) == str)
#     print(type(eg.value) == float)
