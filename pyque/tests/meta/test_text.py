#!/usr/bin/env python3

from unittest import TestCase

from pyque.meta.text import TextStream

instream = """\
&CONTROL
calculation = 'scf'
restart_mode = 'from_scratch'
tstress = .true.
tprnfor = .true.
prefix = 'Fe-hcp'
pseudo_dir = './pseudo'
outdir = './tmp'
etot_conv_thr = 1.0E-6
forc_conv_thr = 1.0D-4
dt = 15
nstep = 100
/
&SYSTEM
ibrav = 0
celldm(1) = 4.493745
nat = 2
ntyp = 1
ecutwfc = 80.0
ecutrho = 320.0
occupations = 'smearing'
degauss = 0.005
smearing = 'mp'
/
&ELECTRONS
diagonalization = 'david'
mixing_mode = 'plain'
mixing_beta = 0.3
conv_thr = 1.0d-8
/
CELL_PARAMETERS
          1.0000000000         0.0000000000         0.0000000000 
         -0.5000000000         0.8660254000         0.0000000000 
          0.0000000000         0.0000000000         1.6040000000  
ATOMIC_SPECIES
Fe   1.0   Fe.KS.GGA-PBE-paw.UPF   
ATOMIC_POSITIONS { crystal }
Fe   0.333333333333333   0.666666666666667   0.250000000000000   
Fe   0.666666666666667   0.333333333333333   0.750000000000000   
K_POINTS { automatic }
4 4 4 1 1 1
"""


class StreamTester(TestCase):
    def setUp(self):
        self.ts = TextStream(instream=instream)

    def test_contents(self):
        self.assertEqual(self.ts.contents, instream)
