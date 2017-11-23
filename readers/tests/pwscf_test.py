#!/usr/bin/env python3
# created at Nov 17, 2017 5:43 PM by Qi Zhang

import unittest

from readers.pwscf import *
from readers.vcrelax import *


class TestPWscfInputReader(unittest.TestCase):
    def setUp(self):
        self.pir = SCFInputReader('test')

    def test_build_input_object(self):
        print('The object is:\n{0}'.format(self.pir.build_input_tree()))


# class TestVCRelaxInputReader(unittest.TestCase):
#     def setUp(self):
#         self.vir = VCRelaxInputfileReader('Fe.in')
#
#     def test_build_vc_relax_input_tree(self):
#         print(self.vir.build_input_tree())
