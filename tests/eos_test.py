#!/usr/bin/env python3
# created at Aug 17, 2017 1:02 PM by Nil-Zil
"""
This is a unit test module of eos.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

from calculators.eos import EOS


class TestEOS(unittest.TestCase):
    def setUp(self):
        self.eos = EOS()

    def test_vinet(self):
        self.assertAlmostEqual(self.eos.vinet(
            100, 137.6852, 283.29, 4.86), 190.99334615)

    def test_solve_vinet(self):
        self.assertAlmostEqual(self.eos.solve_vinet(0, 137.6852, 283.29, 4.86)[0], 137.6852)


if __name__ == "__main__":
    unittest.main()
