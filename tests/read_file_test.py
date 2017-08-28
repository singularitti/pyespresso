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

    def test_read_phonon_dispersion(self):
        self.rpo.read_phonon_dispersion()
