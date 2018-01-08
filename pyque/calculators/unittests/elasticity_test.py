#!/usr/bin/env python3
# created at Oct 20, 2017 11:15 PM by Qi Zhang
"""
This is a unit test module of elasticity.py.
If all tests are OK, it does not mean that the program has no bug.
But if unittest does not pass, there must be bug(s) in the program.
"""

from calculators.elasticity import *

import unittest


class ElasticityCalculationTester(unittest.TestCase):
    def setUp(self):
        self.ec = ElasticityCalculator('../../readers/tests/Cij_file.dat')

    def test_create_compliance_tensor(self):
        print(self.ec.compliance_tensors[0])

    def test_derive_bulk_modulus_voigt_average(self):
        print(self.ec.derive_bulk_modulus_voigt_average(self.ec.elastic_tensors[0]))

    def test_derive_bulk_modulus_reuss_average(self):
        print(self.ec.derive_bulk_modulus_reuss_average(self.ec.compliance_tensors[0]))

    def test_derive_shear_modulus_voigt_average(self):
        print(self.ec.derive_shear_modulus_voigt_average(self.ec.elastic_tensors[0]))

    def test_derive_shear_modulus_reuss_average(self):
        print(self.ec.derive_shear_modulus_reuss_average(self.ec.compliance_tensors[0]))

    def test_derive_bulk_modulus_vrh_average(self):
        print(self.ec.derive_bulk_modulus_vrh_average(self.ec.elastic_tensors[0], self.ec.compliance_tensors[0]))

    def test_derive_shear_modulus_vrh_average(self):
        print(self.ec.derive_shear_modulus_vrh_average(self.ec.elastic_tensors[0], self.ec.compliance_tensors[0]))

    def test_derive_isotropic_poisson_ratio(self):
        print(self.ec.derive_isotropic_poisson_ratio(self.ec.elastic_tensors[0], self.ec.compliance_tensors[0]))

    def test_derive_universal_elastic_anisotropy(self):
        print(self.ec.derive_universal_elastic_anisotropy(self.ec.elastic_tensors[0], self.ec.compliance_tensors[0]))
        self.assertGreaterEqual(
            self.ec.derive_universal_elastic_anisotropy(self.ec.elastic_tensors[0], self.ec.compliance_tensors[0]), 0)


if __name__ == "__main__":
    unittest.main()
