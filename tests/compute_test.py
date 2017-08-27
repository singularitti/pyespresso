#!/usr/bin/env python3
# created at Aug 26, 2017 9:44 PM by Nil-Zil
"""
This is a unit test module of compute.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

import numpy as np

import output.compute as oc


class TestCompute(unittest.TestCase):
    def setUp(self):
        self.cph = oc.ComputePHonon()

    def test_generate_3d_segment(self):
        self.assertTrue(np.array_equal(self.cph.generate_3d_segment([0, 0, 0], [1, 1, 1], 10),
                                       np.repeat(np.linspace(0, 0.9, 10), 3).reshape([10, 3])))
