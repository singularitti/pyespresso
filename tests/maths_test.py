#!/usr/bin/env python3
# created at Aug 26, 2017 9:11 PM by Nil-Zil
"""
This is a unit test module of maths.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

import miscellaneous.maths as mm


class TestMath(unittest.TestCase):
    def test_compute_3d_distance(self):
        print(type(mm.compute_3d_distance([0, 0, 0], [1, 2, 3])))
        self.assertAlmostEqual(mm.compute_3d_distance([0, 0, 0], [1, 2, 3]), 3.7416573868)


if __name__ == "__main__":
    unittest.main()
