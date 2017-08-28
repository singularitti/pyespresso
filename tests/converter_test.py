#!/usr/bin/env python3
# created at Aug 15, 2017 12:24 PM by Nil-Zil
"""
This is a unit test module of converter.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

import unittest

from miscellaneous.converter import *


class TestConverter(unittest.TestCase):
    def test_length(self):
        self.assertAlmostEqual(call_simple_converter(
            'l', 7.4268, 'b', 'a'), 3.930093308203956)

    def test_volume(self):
        self.assertAlmostEqual(call_simple_converter(
            'v', 11.176, 'a3', 'b3'), 75.41938641127665)

    def test_energy(self):
        self.assertAlmostEqual(call_simple_converter(
            'e', 0.02533, 'ry', 'K'), 3999.2920133105276)

    def test_pressure(self):
        self.assertAlmostEqual(call_simple_converter('p', 3, 'mbar', 'gpa'), 300.0)
        self.assertTrue(call_simple_converter('p', [1, 2, 4], 'mbar', 'gpa') == [100.0, 200.0, 400])
        self.assertTrue(np.array_equal(
            call_simple_converter('p', np.array([1, 2, 4]), 'mbar', 'gpa'),
            np.array([100, 200, 400])
        ))


if __name__ == "__main__":
    unittest.main()
