#!/usr/bin/env python3
# created at Nov 20, 2017 3:17 PM by Qi Zhang

from meta.parameter import Parameter, Namelist

INPUTPH_namelist = Namelist(
    'INPUTPH', [
        'amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph', 'alpha_mix', 'nmix_ph', 'verbosity', 'reduce_io',
        'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil', 'lrpa', 'lnoloc', 'trans', 'lraman',
        'eth_rps',
        'eth_ns', 'dek', 'recover', 'low_directory_check', 'only_init', 'qplot', 'q2d', 'q_in_band_form',
        'electron_phonon', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol', 'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir',
        'search_sym', 'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3', 'start_irr', 'last_irr',
        'nat_todo', 'modenum', 'start_q', 'last_q', 'dvscf_star', 'drho_star'
    ], [
        0, './', 'pwscf', 100, 1e-12, 0.7, 4, 'val', False,
        1e7, 'matdyn', '', '', False, False, False, True, False, 1e-9,
        1e-12, 1e-3, False, False, False, False, False, False,
        '', False, False, False, False, False, False, False, False, False, False,
        True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3,
        0, 0, 1, 3, 'disabled', 'disabled'
    ]
)


class INPUTPHParameter(Parameter):
    def __init__(self, name, value):
        default_value, value_type = INPUTPH_namelist.default_parameters[name]
        super().__init__(name, value, value_type)
        self._default_value = default_value
        self._in_namelist = INPUTPH_namelist.__name__
