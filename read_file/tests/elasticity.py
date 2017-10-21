#!/usr/bin/env python3
# created at Oct 20, 2017 11:19 PM by Qi Zhang

import unittest

from read_file.elasticity import *


class TestReadElasticityOutput(unittest.TestCase):
    def setUp(self):
        self.eor = ElasticityOutputReader('Cij_file.dat')

    def test_read_elastic_tensor(self):
        p_list, c_list = self.eor.read_elastic_tensor()
        print(p_list, c_list)


if __name__ == "__main__":
    unittest.main()
