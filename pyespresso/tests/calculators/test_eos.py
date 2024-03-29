#!/usr/bin/env python3
# created at Aug 17, 2017 1:02 PM by Nil-Zil
"""
This is a unit test module of vinet.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

from pyespresso.tools.eos import VinetEoS


class TestVinet(unittest.TestCase):
    def setUp(self):
        self.vinet = VinetEoS()

    def test_vinet(self):
        self.assertAlmostEqual(self.vinet.p_of_v(100, 137.6852, 283.29, 4.86), 190.99334615)

    def test_solve_p_vs_v(self):
        self.assertAlmostEqual(self.vinet.solve_v_by_p(0, 137.6852, 283.29, 4.86)[0], 137.6852)
