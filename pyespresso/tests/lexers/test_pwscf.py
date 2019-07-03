#!/usr/bin/env python3
"""
:mod:`` -- title
========================================

.. module 
   :platform: Unix, Windows, Mac, Linux
   :synopsis: Example is got from
   `here <https://github.com/maxhutch/quantum-espresso/blob/master/PW/examples/example10/run_example>_`.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import unittest

from pyespresso.core.cards import *
from pyespresso.lexers.pwscf import PWscfInputLexer

instream_sample = """\
&control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='silicon',
    lelfield=.true.,
    nberrycyc=1
    pseudo_dir='$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
 /
 &system
    ibrav= 1, celldm(1)=10.18, nat=  8, ntyp= 1,
    ecutwfc = 20.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random',
    efield_cart(1)=0.d0,efield_cart(2)=0.d0,efield_cart(3)=0.d0
 /
ATOMIC_SPECIES
 Si  28.086 Si.pbe-rrkj.UPF
ATOMIC_POSITIONS
 Si -0.125 -0.125 -0.125
 Si  0.375  0.375 -0.125
 Si  0.375 -0.125  0.375
 Si -0.125  0.375  0.375
 Si  0.125  0.125  0.125
 Si  0.625  0.625  0.125
 Si  0.625  0.125  0.625
 Si  0.125  0.625  0.625
K_POINTS {automatic}
3 3 7 0 0 0
"""


class TestPWscfInputLexer(unittest.TestCase):
    def setUp(self):
        self.lexer = PWscfInputLexer(inp=instream_sample)

    def test_namelists_found(self):
        self.assertEqual(self.lexer.namelists_found, {'&CONTROL', '&SYSTEM', '&ELECTRONS'})

    def test_cards_found(self):
        self.assertEqual(self.lexer.cards_found, {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'})

    # def test_get_control_namelist(self):
    #     print(self.__lexer.get_control_namelist())
    #     pass

    def test_lex_control_namelist(self):
        self.assertEqual(self.lexer.lex_control_namelist(),
                         {'calculation': "'scf'", 'restart_mode': "'from_scratch'", 'prefix': "'silicon'",
                          'lelfield': '.true.', 'nberrycyc': '1', 'pseudo_dir': "'$PSEUDO_DIR/'",
                          'outdir': "'$TMP_DIR/'"})

    def test_lex_system_namelist(self):
        self.assertEqual(self.lexer.lex_system_namelist(),
                         {'ibrav': '1', 'celldm(1)': '10.18', 'nat': '8', 'ntyp': '1', 'ecutwfc': '20.0'})

    def test_lex_electrons_namelist(self):
        self.assertEqual(self.lexer.lex_electrons_namelist(),
                         {'diagonalization': "'david'", 'conv_thr': '1.0d-8', 'mixing_beta': '0.5',
                          'startingwfc': "'random'", 'efield_cart(1)': '0.d0', 'efield_cart(2)': '0.d0',
                          'efield_cart(3)': '0.d0'})

    def test_lex_ions_namelist(self):
        self.assertEqual(self.lexer.lex_ions_namelist(), None)

    def test_lex_cell_namelist(self):
        self.assertEqual(self.lexer.lex_cell_namelist(), None)

    def test_lex_atomic_species(self):
        self.assertEqual(self.lexer.lex_atomic_species(),
                         [AtomicSpecies(name='Si', mass='28.086', pseudopotential='Si.pbe-rrkj.UPF')])

    def test_lex_atomic_positions(self):
        self.assertEqual(self.lexer.lex_atomic_positions(),
                         ([AtomicPosition(name='Si', pos=['-0.125', '-0.125', '-0.125'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.375', '0.375', '-0.125'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.375', '-0.125', '0.375'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['-0.125', '0.375', '0.375'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.125', '0.125', '0.125'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.625', '0.625', '0.125'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.625', '0.125', '0.625'], if_pos=['1', '1', '1']),
                           AtomicPosition(name='Si', pos=['0.125', '0.625', '0.625'], if_pos=['1', '1', '1'])],
                          'alat'))

    def test_lex_k_points(self):
        self.assertEqual(self.lexer.lex_k_points(),
                         (AutomaticKPoints(grid=['3', '3', '7'], offsets=['0', '0', '0']), 'automatic'))
