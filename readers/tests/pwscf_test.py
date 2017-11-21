#!/usr/bin/env python3
# created at Nov 17, 2017 5:43 PM by Qi Zhang

import unittest

from readers.pwscf import *
from readers.vcrelax import *


class TestPWscfInputReader(unittest.TestCase):
    def setUp(self):
        self.pir = PWscfInputReader('Fe.in')

    def test_tree(self):
        self.assertTrue(hasattr(self.pir, 'tree'))
        print(self.pir.tree)

    def test_build_PWscf_input_object(self):
        print('The object is:\n{0}'.format(self.pir.build_pwscf_input_object()))

    def test__call__(self):
        print(self.pir())

    def test__str__(self):
        print(self.pir)


class TestVCRelaxInputReader(unittest.TestCase):
    def setUp(self):
        self.vir = VCRelaxInputfileReader('Fe.in')

    def test_build_vc_relax_input_tree(self):
        print(self.vir.build_vc_relax_input_tree())
