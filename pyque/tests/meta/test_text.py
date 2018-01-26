#!/usr/bin/env python3
"""
The example of Quantum ESPRESSO input is got from
`here <https://github.com/maxhutch/quantum-espresso/blob/master/PW/examples/example02/run_example>`_.
"""

import unittest

from pyque.meta.text import TextStream

inp = """\
&CONTROL
  calculation  = "relax",
  prefix       = "CO",
  pseudo_dir   = "$PSEUDO_DIR",
  outdir       = "$TMP_DIR",
/
&SYSTEM
  ibrav     = 0,
  nat       = 2,
  ntyp      = 2,
  ecutwfc   = 24.D0,
  ecutrho   = 144.D0,
/
&ELECTRONS
  conv_thr    = 1.D-7,
  mixing_beta = 0.7D0,
/
&IONS
/
CELL_PARAMETERS bohr
12.0  0.0  0.0
 0.0 12.0  0.0
 0.0  0.0 12.0
ATOMIC_SPECIES
O  1.00  O.pz-rrkjus.UPF
C  1.00  C.pz-rrkjus.UPF
ATOMIC_POSITIONS {bohr}
C  2.256  0.0  0.0
O  0.000  0.0  0.0  0 0 0
K_POINTS {Gamma}
"""


class StreamTester(unittest.TestCase):
    def setUp(self):
        self.ts = TextStream(inp)

    def test_contents(self):
        self.assertEqual(self.ts.contents, inp)
