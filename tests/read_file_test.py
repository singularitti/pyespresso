#!/usr/bin/env python3
# created at Aug 27, 2017 10:15 PM by Nil-Zil
"""
This is a unit test module of read_file.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

from read_file.read_file import *


class TestReadPHononOutput(unittest.TestCase):
    def setUp(self):
        self.rpw = ReadPWscfOutput()
        self.rph = ReadPHononOutput()
        # self.q_dict = self.test_read_q_points()

        # def test_read_k_mesh(self):
        #     print(self.rpw.read_k_mesh('pw_data/eg0.in'))

        # def test_read_q_points(self):
        #     return self.rph.read_q_points('phonon_data/high_sym_q')
        #
        # def test_read_phonon_dispersion(self):
        #     self.rph.read_phonon_dispersion('phonon_data/freq.out', 'Γ->M->K->Γ->A->K')


class TestReadElasticityOutput(unittest.TestCase):
    def setUp(self):
        self.reo = ReadElasticityOutput()

    def test_read_elastic_tensor(self):
        p_list, c_list = self.reo.read_elastic_tensor('cij_data/Cij_file.dat')
        print(p_list, c_list)


if __name__ == "__main__":
    unittest.main()
