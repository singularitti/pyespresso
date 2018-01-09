#!/usr/bin/env python3
# created at Oct 19, 2017 1:03 AM by Qi Zhang
"""
Parameters of pw.x that controls the whole calculation.
Here I construct different sets of parameters, for consideration of performance.
Since finding an element in a list is O(N) but in a set is O(1).
"""

from pyque.meta.parameter import ParameterGeneric, Namelist

# =================================== I am a cut line ===================================
CONTROL_NAMELIST: Namelist = Namelist(
    'CONTROL', [
        'calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress', 'tprnfor',
        'dt', 'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr', 'disk_io',
        'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm', 'lberry', 'gdir', 'nppstr', 'lfcpopt',
        'gate'
    ], [
        'scf', ' ', 'low', 'from_scratch', True, 1, 1, False, False,
        20.0e0, './', './', 'pwscf', True, 1.0e7, 1.0e-4, 1.0e-3, 'medium',
        '$ESPRESSO_PSEUDO', False, False, False, 1, False, False, 1, 1, False, False
    ]
)

# =================================== I am a cut line ===================================
SYSTEM_NAMELIST: Namelist = Namelist(
    'SYSTEM', [
        'ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'tot_charge',
        'starting_charge',
        'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3', 'nr1s',
        'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac', 'occupations',
        'one_atom_occupations', 'starting_spin_angle', 'degauss', 'smearing', 'nspin', 'noncolin', 'ecfixed', 'qcutz',
        'q2sigma', 'input_dft', 'exx_fraction', 'screening_parameter', 'exxdiv_treatment', 'x_gamma_extrapolation',
        'ecutvcut', 'nqx1', 'nqx2', 'nqx3', 'lda_plus_u', 'lda_plus_u_kind', 'Hubbard_U', 'Hubbard_J0', 'Hubbard_alpha',
        'Hubbard_beta', 'Hubbard_J', 'starting_ns_eigenvalue', 'U_projection_type', 'edir', 'emaxpos',
        'eopreg', 'eamp', 'angle1', 'angle2', 'constrained_magnetization', 'fixed_magnetization', 'lambda', 'report',
        'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit', 'fcp_mu', 'vdw_corr', 'london',
        'london_s6', 'london_c6', 'london_rvdw', 'london_rcut', 'ts_vdw_econv_thr', 'ts_vdw_isolated', 'xdm', 'xdm_a1',
        'xdm_a2', 'space_group', 'uniqueb', 'origin_choice', 'rhombohedral', 'zgate', 'relaxz', 'block', 'block_1',
        'block_2', 'block_height'
    ], [
        0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 1, 1, 20, 0.0,
        0.0,
        -1.0, 0.0, 90.0, 360.0, 120.0, 24, 24, 24, 24,
        24, 24, False, False, False, False, False, False, 'smearing',
        False, False, 0.0e0, 'gaussian', 1, False, 0.0, 0.0,
        0.1, 'PBE0', 0.25, 0.106, 'gygi-baldereschi', True,
        0.0, 1, 1, 1, False, 0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, -1.0e0, 'atomic', 1, 0.5e0,
        0.1e0, 0.001, 0.0, 0.0, 'none', 0.0e0, 1.0e0, 100,
        False, 'none', 'pbc', 0.0e0, 0.0e0, 4, 0.0e0, 'none', False,
        0.75, 0.0e0, 0.0e0, 200, 1.0e-6, False, False, 0.6836,
        1.5045, 0, False, 1, True, 0.5, False, False, 0.45,
        0.55, 0.1
    ]
)

# =================================== I am a cut line ===================================
ELECTRONS_NAMELIST: Namelist = Namelist(
    'ELECTRONS', [
        'electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init', 'conv_thr_multi',
        'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns', 'diagonalization', 'ortho_para',
        'diago_thr_init',
        'diago_cg_maxiter', 'diago_david_ndim', 'diago_full_acc', 'efield', 'efield_cart', 'efield_phase',
        'startingpot', 'startingwfc', 'tqr'
    ], [
        100, True, 1.0e-6, False, 1.0e-3, 1.0e-1,
        'plain', 0.7e0, 8, 0, 'david', 0,
        1.0e-6,
        400, 4, False, 0.0e0, (0.0e0, 0.0e0, 0.0e0), 'none',
        'atomic', 'atomic+random', False
    ]
)

# =================================== I am a cut line ===================================
IONS_NAMELIST: Namelist = Namelist(
    'IONS', [
        'ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot',
        'ion_temperature',
        'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim', 'trust_radius_max',
        'trust_radius_min', 'trust_radius_ini', 'w_1', 'w_2'
    ], [
        'bfgs', 'default', 'atomic', 'none', False,
        'not_controlled',
        300.0e0, 100.0e0, 1.0e0, 1, False, 100.0e0, 1, 0.8e0,
        1.0e-3, 0.5e0, 0.01e0, 0.5e0
    ]
)

# =================================== I am a cut line ===================================
CELL_NAMELIST: Namelist = Namelist(
    'CELL', [
        'cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree'
    ], [
        'none', 0.0e0, 0.001, 2.0, 0.5e0, 'all'
    ]
)

# =================================== I am a cut line ===================================
INPUTPH_NAMELIST = Namelist(
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


# =================================== I am a cut line ===================================
class _Parameter(ParameterGeneric):
    """
    This is a prototype for a parameter.
    You only need to provide the name of your parameter, the value of it, and the parameter ``dict`` it belongs to.
    """

    def __init__(self, name: str, value: str, namelist: Namelist):
        self._default_value, value_type = namelist.default_parameters[name]
        super().__init__(name, value, value_type)
        self._in_namelist = namelist.__name__


class CONTROLParameter(_Parameter):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CONTROL_NAMELIST)


class SYSTEMParameter(_Parameter):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, SYSTEM_NAMELIST)


class ELECTRONSParameter(_Parameter):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, ELECTRONS_NAMELIST)


class IONSParameter(_Parameter):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, IONS_NAMELIST)


class CELLParameter(_Parameter):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CELL_NAMELIST)


class INPUTPHParameter(_Parameter):
    def __init__(self, name, value):
        super().__init__(name, value, INPUTPH_NAMELIST)
