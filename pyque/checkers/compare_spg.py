#!/usr/bin/env python3
"""
:mod:`` -- 
========================================

.. module 
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from typing import *
import pyque.core.cell as _cell
from pyque.core.qe_input import VCRelaxInput


def compare_spacegroup(vc_relax_input: VCRelaxInput, vc_relax_output: object) -> bool:
    lattice = vc_relax_input.cell_parameters
    positions = vc_relax_input.atomic_positions
    numbers = vc_relax_input.system_namelist.ntyp
    sym_ops: int = vc_relax_output.get_symmetry_opreations_number()
    cell = _cell.Cell(lattice=lattice, positions=positions, numbers=numbers)
    if len(cell.rotations) == sym_ops:
        return True
    else:
        return False
