#!/usr/bin/env python3
# created at Oct 19, 2017 1:03 AM by Qi Zhang
"""
Parameters of pw.x that controls the whole calculation.
Here I construct different sets of parameters, for consideration of performance.
Since finding an element in a list is O(N) but in a set is O(1).
"""

from miscellaneous.dictionary import merge_dicts

CONTROL = {'calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress', 'tprnfor',
           'dt', 'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr',
           'disk_io', 'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm', 'lberry', 'gdir', 'nppstr',
           'lfcpopt', 'monopole'}

SYSTEM = {'ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'tot_charge',
          'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3', 'nr1s',
          'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac', 'occupations',
          'one_atom_occupations', 'starting_spin_angle', 'degauss', 'smearing', 'nspin', 'noncolin', 'ecfixed', 'qcutz',
          'q2sigma', 'input_dft', 'exx_fraction', 'screening_parameter', 'exxdiv_treatment', 'x_gamma_extrapolation',
          'ecutvcut', 'nqx1', 'nqx2', 'nqx3', 'lda_plus_u', 'lda_plus_u_kind', 'Hubbard_U', 'Hubbard_J0',
          'Hubbard_alpha', 'Hubbard_beta', 'Hubbard_J(i,ityp)', 'starting_ns_eigenvalue(m,ispin,I)',
          'U_projection_type', 'edir', 'emaxpos', 'eopreg', 'eamp', 'angle1', 'angle2', 'constrained_magnetization',
          'fixed_magnetization', 'lambda', 'report', 'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield',
          'esm_nfit', 'fcp_mu', 'vdw_corr', 'london', 'london_s6', 'london_c6', 'london_rvdw', 'london_rcut',
          'ts_vdw_econv_thr', 'ts_vdw_isolated', 'xdm', 'xdm_a1', 'xdm_a2', 'space_group', 'uniqueb', 'origin_choice',
          'rhombohedral', 'zmon', 'realxz', 'block', 'block_1', 'block_2', 'block_height'}

ELECTRONS = {'electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init', 'conv_thr_multi',
             'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns', 'diagonalization', 'ortho_para',
             'diago_thr_init', 'diago_cg_maxiter', 'diago_david_ndim', 'diago_full_acc', 'efield', 'efield_cart',
             'efield_phase', 'startingpot', 'startingwfc', 'tqr'}

IONS = {'ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot',
        'ion_temperature', 'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim',
        'trust_radius_max', 'trust_radius_min', 'trust_radius_ini', 'w_1', 'w_2'}

CELL = {'cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree'}

default_control = {'restart_mode': 'from_scratch',
                   'tstress': '.true.',
                   'tprnfor': '.true.',
                   'verbosity': 'high'}

default_system = {'ibrav': '0'}

default_electrons = {'diagonalization': 'david'}

default_parameters = merge_dicts(default_control, default_system, default_electrons)
