#!/usr/bin/env python3
"""
:mod:`mod` -- title
========================================

.. module mod
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from typing import *

supported_content: Dict[str, str] = {
    'lattice_parameter': "lattice parameter \(alat\)\s+=\s*(\d+\.\d+)",
    'total_energy': "!\s+total\s+energy\s+=\s+(-?\d+\.\d+)",
    'cell_volume': "unit-cell volume\s+=\s+(\d+\.\d+)",
    'pressure': "P=\s+(-?\d+\.\d+)",
    'kinetic_energy_cutoff': "kinetic-energy cutoff\s+=\s*(-?\d+\.\d+)",
    'charge_density_cutoff': "charge density cutoff\s+=\s*(-?\d+\.\d+)",
    'atoms_num_per_cell': "number of atoms\/cell\s+=\s*(\d+)",
    'atoms_types_num': "number of atomic INPUTPH_types\s+=\s*(\d+)",
    'electrons_num': "number of electrons\s+=\s*(-?\d+\.\d+)",
    'ks_states_num': "number of Kohn-Sham states\s*=\s*(\d+)",
    'mixing_beta': "mixing beta\s+=\s*(-?\d*\.\d+)",
    'nstep': "nstep\s+=\s*(\d+)",
    'iteration_num': "number of iterations used\s+=\s*(\d+)",
    'symmetry_operations_num': "(\d+)\s+Sym\. Ops.*found",
}

stress_mapping: Dict[str, str] = {
    'kinetic': 'kinetic',
    'local': 'local',
    'nonlocal': 'nonloc.',
    'hartree': 'hartree',
    'exc-cor': 'exc-cor',
    'corecor': 'corecor',
    'ewald': 'ewald',
    'hubbard': 'hubbard',
    'london': 'london',
    'XDM': 'XDM',
    'dft-nl': 'dft-nl',
    'TS-vdW': 'TS-vdW'
}
