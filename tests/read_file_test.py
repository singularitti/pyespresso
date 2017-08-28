#!/usr/bin/env python3
# created at Aug 27, 2017 10:15 PM by Nil-Zil
"""
This is a unit test module of read_file.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

from output.read_file import *


class TestReadPHononOutput(unittest.TestCase):
    def setUp(self):
        self.rpo = ReadPHononOutput()
        self.q_dict = self.test_read_q_points()

    def test_read_q_points(self):
        return self.rpo.read_q_points('phonon_data/high_sym_q')

    def test_read_phonon_dispersion(self):
        self.rpo.read_phonon_dispersion('phonon_data/freq.out', self.q_dict, 'GM->M->K->GM->A->K')


if __name__ == "__main__":
    unittest.main()
