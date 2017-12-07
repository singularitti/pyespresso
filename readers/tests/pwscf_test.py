#!/usr/bin/env python3
# created at Nov 17, 2017 5:43 PM by Qi Zhang

import unittest

from basics.pw_builder import *


class TestPWscfInputReader(unittest.TestCase):
    def setUp(self):
        self.pwsi = build_pw_input('test')

    def test_write_to_file(self):
        self.pwsi.write_to_file('new')
