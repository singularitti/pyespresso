#!/usr/bin/env python3
# created at Nov 25, 2017 8:05 PM by Qi Zhang
import unittest

from readers.phonon import *


class TestPHononInputReader(unittest.TestCase):
    def setUp(self):
        self.phir = PhononInputReader('ph.in')

    def test_build_phono_input_tree(self):
        print(self.phir.build_input_tree())
