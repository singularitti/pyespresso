#!/usr/bin/env python3
# created at Aug 26, 2017 9:44 PM by Nil-Zil
"""
This is a unit test module of compute.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

import numpy as np

import miscellaneous.compute as oc


class TestCompute(unittest.TestCase):
    def setUp(self):
        self.cph = oc.ComputePHonon()

    def test_generate_3d_segment(self):
        m = np.array([[2, 0, 1],
                      [1.9, 0, 1],
                      [1.8, 0, 1],
                      [1.7, 0, 1],
                      [1.6, 0, 1],
                      [1.5, 0, 1],
                      [1.4, 0, 1],
                      [1.3, 0, 1],
                      [1.2, 0, 1],
                      [1.1, 0, 1],
                      [1, 0, 1]])
        self.assertTrue(np.allclose(self.cph.generate_3d_segment([2, 0, 1], [1, 0, 1], 11), m))

    def test_generate_q_path(self):
        self.cph.generate_q_path('phonon_data/high_sym_q', 'phonon_data/qpts', 'Γ->M->K->Γ->A->K', 100)


if __name__ == "__main__":
    unittest.main()
