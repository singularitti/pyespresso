#!/usr/bin/env python3
# created at Nov 17, 2017 5:43 PM by Qi Zhang

import unittest

from beeprint import pp

from basics.pw_builder import *
from basics.pwscf import *


# pp = pprint.PrettyPrinter(indent=4)


class TestPWscfInputReader(unittest.TestCase):
    def setUp(self):
        self.pwsi = build_pw_input('test')
        # pp(self.pwsi)
        self.fancy = PWInputFancier(self.pwsi)

    def test_write_to_file(self):
        self.pwsi.write_to_file('new')

    def test_print_pw_input(self):
        print_pw_input(self.pwsi)

    def test_fancy_CONTROL(self):
        s = self.fancy.fancy_CONTROL()
        pp(s)

    def test_fancy_ELECTRONS(self):
        s = self.fancy.fancy_ELECTRONS()
        print(s)

    def test_fancy_input(self):
        pp(self.fancy.fancy_input())
