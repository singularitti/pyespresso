#!/usr/bin/env python3
# created at Oct 19, 2017 1:03 AM by Qi Zhang
"""
Parameters of pw.x that controls the whole calculation.
Here I construct different sets of parameters, for consideration of performance.
Since finding an element in a list is O(N) but in a set is O(1).
"""

from collections import defaultdict
from typing import DefaultDict, Type, Union, Tuple

from basics.parameter import Parameter, Namelist

# Type alias
DefaultParameters = DefaultDict[str, Tuple[Union[str, int, float, bool], Type[Union[str, int, float, bool]]]]

# This is a list of names for 'CONTROL' namelist.
CONTROL_names = [
    'calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress', 'tprnfor', 'dt',
    'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr', 'disk_io',
    'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm', 'lberry', 'gdir', 'nppstr', 'lfcpopt', 'gate'
]

# This is a list of QE default values for 'CONTROL' namelist.
CONTROL_default_values = [
    'scf', ' ', 'low', 'from_scratch', True, 1, 1, False, False, 20.0e0,
    './', './', 'pwscf', True, 1.0e7, 1.0e-4, 1.0e-3, 'medium',
    '$ESPRESSO_PSEUDO', False, False, False, 1, False, False, 1, 1, False, False
]

# This is a list of types of each default value for 'CONTROL' namelist.
CONTROL_value_types = [type(x) for x in CONTROL_default_values]

# This is a `DefaultDict` of `(QE default value, types of each default value)` tuples for 'CONTROL' namelist, with those
# names to be its keys.
CONTROL_default_parameters: DefaultParameters = defaultdict(tuple)
for i, k in enumerate(CONTROL_names):
    CONTROL_default_parameters[k] = (CONTROL_default_values[i], CONTROL_value_types[i])

# This is a `NamedTuple` for 'CONTROL' namelist.
CONTROLNamelist = Namelist(name='CONTROL', keys=set(CONTROL_names))
# =================================== I am a cut line ===================================

SYSTEM_names = [
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

SYSTEM_default_values = [0 for _ in range(len(SYSTEM_names))]

SYSTEM_default_parameters: DefaultParameters = dict(zip(SYSTEM_names, SYSTEM_default_values))
for i, k in enumerate(SYSTEM_default_parameters):
    SYSTEM_default_parameters[k] = (SYSTEM_default_values[i], SYSTEM_default_parameters[i])

SYSTEMNamelist = Namelist(name='SYSTEM', keys=set(SYSTEM_names))
# =================================== I am a cut line ===================================

ELECTRONS_names = [
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

ELECTRONS_default_parameters: DefaultParameters = dict(zip(ELECTRONS_names, ELECTRONS_default_values))
for i, k in enumerate(ELECTRONS_names):
    ELECTRONS_default_parameters[k] = (ELECTRONS_default_values[i], ELECTRONS_value_types[i])

ELECTRONSNamelist = Namelist(name='ELECTRONS', keys=set(ELECTRONS_names))
# =================================== I am a cut line ===================================

IONS_names = [
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

IONS_default_parameters: DefaultParameters = defaultdict(tuple)
for i, k in enumerate(IONS_names):
    IONS_default_parameters[k] = (IONS_default_values[i], IONS_value_types[i])

IONSNamelist = Namelist(name='IONS', keys=set(IONS_names))
# =================================== I am a cut line ===================================

CELL_names = ['cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree']

CELL_default_values = ['none', 0.0e0, 0.001, 2.0, 0.5e0, 'all']

CELL_value_types = [type(x) for x in CELL_default_values]

CELL_default_parameters: DefaultParameters = defaultdict(tuple)
for i, k in enumerate(CELL_names):
    CELL_default_parameters[k] = (CELL_default_values[i], CELL_value_types[i])

CELLNamelist = Namelist(name='CELL', keys=set(CELL_names))


# =================================== I am a cut line ===================================
class PWParameterGeneric(Parameter):
    """
    This is a generic for building a parameter for one of 'CONTROL', 'SYSTEM', `ELECTRONS`, `IONS`, or `CELL` namelists.
    You only need to provide the name of your parameter, the value of it, and the parameter `dict` it belongs to.
    """

    def __init__(self, name: str, value: str, default_parameters: dict):
        default_parameter = default_parameters[name]
        super().__init__(name, value, default_parameter[1])
        self._default_value = default_parameter[0]

    @property
    def default_value(self):
        """
        Return the default value defined by Quantum ESPRESSO.

        :return: The default value defined by Quantum ESPRESSO.
        """
        return self._default_value


class CONTROLParameter(PWParameterGeneric):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CONTROL_default_parameters)


class SYSTEMParameter(PWParameterGeneric):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, SYSTEM_default_parameters)


class ELECTRONSParameter(PWParameterGeneric):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, ELECTRONS_default_parameters)


class IONSParameter(PWParameterGeneric):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, IONS_default_parameters)


class CELLParameter(PWParameterGeneric):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CELL_default_parameters)
