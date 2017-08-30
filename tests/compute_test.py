#!/usr/bin/env python3
# created at Aug 26, 2017 9:44 PM by Nil-Zil
"""
This is a unit test module of compute.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

import numpy as np

import miscellaneous.compute as mc


class TestCompute(unittest.TestCase):
    def setUp(self):
        self.cph = mc.ComputePHonon()
        self.gp = mc.ReciprocalPathGenerator('phonon_data/high_sym_q', 'Γ->M->K->Γ->A->K')

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
        self.assertTrue(np.allclose(self.gp.linspace_3d([2, 0, 1], [1, 0, 1], 11), m))

    def test_generate_q_path(self):
        # self.gp.generate_q_path(100)
        print(self.gp.generate_q_path([4, 4, 4, 3, 3]))
        # self.gp.generate_q_path(np.array([100, 100, 200, 50, 100]))


if __name__ == "__main__":
    unittest.main()
