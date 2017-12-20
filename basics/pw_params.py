#!/usr/bin/env python3
# created at Oct 19, 2017 1:03 AM by Qi Zhang
"""
Parameters of pw.x that controls the whole calculation.
Here I construct different sets of parameters, for consideration of performance.
Since finding an element in a list is O(N) but in a set is O(1).
"""

from basics.parameter import Parameter, Namelist

CONTROL_names = [
    'calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress', 'tprnfor', 'dt',
    'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr', 'disk_io',
    'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm', 'lberry', 'gdir', 'nppstr', 'lfcpopt', 'gate'
]

CONTROL_default_values = [
    'scf', ' ', 'low', 'from_scratch', True, 1, 1, False, False, 20.0e0,
    './', './', 'pwscf', True, 1.0e7, 1.0e-4, 1.0e-3, 'medium',
    '$ESPRESSO_PSEUDO', False, False, False, 1, False, False, 1, 1, False, False
]

CONTROL_value_types = [type(x) for x in CONTROL_default_values]

_CONTROL_default_parameters = dict(zip(CONTROL_names, CONTROL_default_values))

CONTROLNamelist = Namelist(name='CONTROL', keys=CONTROL_names)
# =================================== I am a cut line ===================================

SYSTEM_keys = [
    'ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'tot_charge', 'starting_charge',
    'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3', 'nr1s',
    'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac', 'occupations',
    'one_atom_occupations', 'starting_spin_angle', 'degauss', 'smearing', 'nspin', 'noncolin', 'ecfixed', 'qcutz',
    'q2sigma', 'input_dft', 'exx_fraction', 'screening_parameter', 'exxdiv_treatment', 'x_gamma_extrapolation',
    'ecutvcut', 'nqx1', 'nqx2', 'nqx3', 'lda_plus_u', 'lda_plus_u_kind', 'Hubbard_U', 'Hubbard_J0', 'Hubbard_alpha',
    'Hubbard_beta', 'Hubbard_J(i,ityp)', 'starting_ns_eigenvalue(m,ispin,I)', 'U_projection_type', 'edir', 'emaxpos',
    'eopreg', 'eamp', 'angle1', 'angle2', 'constrained_magnetization', 'fixed_magnetization', 'lambda', 'report',
    'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit', 'fcp_mu', 'vdw_corr', 'london',
    'london_s6', 'london_c6', 'london_rvdw', 'london_rcut', 'ts_vdw_econv_thr', 'ts_vdw_isolated', 'xdm', 'xdm_a1',
    'xdm_a2', 'space_group', 'uniqueb', 'origin_choice', 'rhombohedral', 'zgate', 'relaxz', 'block', 'block_1',
    'block_2', 'block_height'
]

SYSTEM_default_values = [0 for _ in range(len(SYSTEM_keys))]

_SYSTEM_default_parameters = dict(zip(SYSTEM_keys, SYSTEM_default_values))

SYSTEMNamelist = Namelist(name='SYSTEM', keys=SYSTEM_keys)
# =================================== I am a cut line ===================================

ELECTRONS_keys = [
    'electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init', 'conv_thr_multi',
    'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns', 'diagonalization', 'ortho_para', 'diago_thr_init',
    'diago_cg_maxiter', 'diago_david_ndim', 'diago_full_acc', 'efield', 'efield_cart', 'efield_phase', 'startingpot',
    'startingwfc', 'tqr'
]

ELECTRONS_default_values = [
    100, True, 1.0e-6, False, 1.0e-3, 1.0e-1,
    'plain', 0.7e0, 8, 0, 'david', 0, 1.0e-6,
    400, 4, False, 0.0e0, (0.0e0, 0.0e0, 0.0e0), 'none', 'atomic',
    'atomic+random', False
]

ELECTRONS_value_types = [type(x) for x in ELECTRONS_default_values]

_ELECTRONS_default_parameters = dict(zip(ELECTRONS_keys, ELECTRONS_default_values))

ELECTRONSNamelist = Namelist(name='ELECTRONS', keys=ELECTRONS_keys)
# =================================== I am a cut line ===================================

IONS_keys = [
    'ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot', 'ion_temperature',
    'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim', 'trust_radius_max', 'trust_radius_min',
    'trust_radius_ini', 'w_1', 'w_2'
]

IONS_default_values = [
    'bfgs', 'default', 'atomic', 'none', False, 'not_controlled',
    300.0e0, 100.0e0, 1.0e0, 1, False, 100.0e0, 1, 0.8e0, 1.0e-3,
    0.5e0, 0.01e0, 0.5e0
]

IONS_value_types = [type(x) for x in IONS_default_values]

_IONS_default_parameters = dict(zip(IONS_keys, IONS_default_values))

IONSNamelist = Namelist(name='IONS', keys=IONS_keys)
# =================================== I am a cut line ===================================

CELL_keys = ['cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree']

CELL_default_values = ['none', 0.0e0, 0.001, 2.0, 0.5e0, 'all']

CELL_value_types = [type(x) for x in CELL_default_values]

_CELL_default_parameters = dict(zip(CELL_keys, CELL_default_values))

CELLNamelist = Namelist(name='CELL', keys=CELL_keys)


class PWParameterGeneric(Parameter):
    def __init__(self, name: str, value: str, default_parameters: dict):
        super().__init__(name, value)
        self._default_value = default_parameters[name]

    @property
    def default_value(self):
        return self._default_value


class CONTROLParameter(PWParameterGeneric):
    def __init__(self, name: str, value: str):
        super().__init__(name, value, _CONTROL_default_parameters)


class SYSTEMParameter(PWParameterGeneric):
    def __init__(self, name: str, value: str):
        super().__init__(name, value, _SYSTEM_default_parameters)


class ELECTRONSParameter(PWParameterGeneric):
    def __init__(self, name: str, value: str):
        super().__init__(name, value, _ELECTRONS_default_parameters)


class IONSParameter(PWParameterGeneric):
    def __init__(self, name: str, value: str):
        super().__init__(name, value, _IONS_default_parameters)


class CELLParameter(PWParameterGeneric):
    def __init__(self, name: str, value: str):
        super().__init__(name, value, _CELL_default_parameters)
