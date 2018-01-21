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

from pyque.lexer.pwscf import PWscfInputLexer

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
        self.__lexer = PWscfInputLexer(instream=instream_sample)

    def test_namelists_found(self):
        self.assertEqual(self.__lexer.namelists_found, {'&CONTROL', '&SYSTEM', '&ELECTRONS'})

    def test_cards_found(self):
        self.assertEqual(self.__lexer.cards_found, {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'})

    def test_get_control_namelist(self):
        # print(self.__lexer.get_control_namelist())
        pass

    def test_lex_control_namelist(self):
        print(repr(self.__lexer.lex_control_namelist()))
        # self.assertEqual(self.__lexer.lex_control_namelist(),
        #                  {'calculation': "calculation 'scf' str scf CONTROL",
        #                   'restart_mode': "restart_mode 'from_scratch' str from_scratch CONTROL",
        #                   'prefix': "prefix 'silicon' str pwscf CONTROL",
        #                   'lelfield': "lelfield.true.bool False CONTROL", 'nberrycyc': "nberrycyc 1 int 1 CONTROL",
        #                   'pseudo_dir': "pseudo_dir '$PSEUDO_DIR/' str $ESPRESSO_PSEUDO CONTROL",
        #                   'outdir': "outdir '$TMP_DIR/' str. / CONTROL"})

    def test_lex_atomic_species(self):
        print(self.__lexer.lex_atomic_species())
        self.assertEqual(self.__lexer.lex_atomic_species(),
                         [AtomicSpecies(name='Si', mass='28.086', pseudopotential='Si.pbe-rrkj.UPF')])
