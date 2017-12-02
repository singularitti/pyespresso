#!/usr/bin/env python3
# created at Dec 2, 2017 2:58 AM by Qi Zhang

import unittest

from calculators.spacegroup import *
from calculators.unit_tests.spacegroup_test_data import *

rutile = Cell(*rutile_spglib_cell)


class CellTester(unittest.TestCase):
    def test_get_spacegroup(self):
        self.assertEqual(rutile.get_spacegroup(), 'P4_2/mnm (136)')

    def test_get_symmetry(self):
        rs = rutile.get_symmetry()
        for k in rutile_symmetry.keys():
            np.testing.assert_array_equal(rs[k], rutile_symmetry[k])

    # def test_standardize_cell(self):

    def test_spacegroup(self):
        self.assertEqual(rutile.spacegroup, 'P4_2/mnm (136)')

    def test_symmetry(self):
        rs = rutile.symmetry
        for k in rutile_symmetry.keys():
            np.testing.assert_array_equal(rs[k], rutile_symmetry[k])

    def test_number(self):
        self.assertEqual(rutile.number, 136)

    def test_hall_number(self):
        self.assertEqual(rutile.hall_number, 419)

    def test_international(self):
        self.assertEqual(rutile.international, 'P4_2/mnm')

    def test_hall(self):
        self.assertEqual(rutile.hall, '-P 4n 2n')

    def test_choice(self):
        self.assertEqual(rutile, '')

    def test_transformation_matrix(self):
        np.testing.assert_array_equal(rutile.transformation_matrix, rutile_transformation_matrix)

    def test_origin_shift(self):
        np.testing.assert_array_equal(rutile.origin_shift, np.array([0., 0., 0.]))

    def test_rotations(self):
        np.testing.assert_array_equal(rutile.rotations, rutile_symmetry['rotations'])

    def test_translations(self):
        np.testing.assert_array_equal(rutile.translations, rutile_symmetry['translations'])

    def test_wyckoffs(self):
        self.assertEqual(rutile.wyckoffs, ['a', 'a', 'f', 'f', 'f', 'f'])

    def test_equivalent_atoms(self):
        np.testing.assert_array_equal(rutile.equivalent_atoms, rutile_symmetry['equivalent_atoms'])

    def test_mapping_to_primitive(self):
        np.testing.assert_array_equal(rutile.mapping_to_primitive, np.array([0, 1, 2, 3, 4, 5], dtype=np.int32))

    def test_std_lattice(self):
        np.testing.assert_array_equal(rutile.std_lattice, rutile_dataset['std_lattice'])

    def test_std_types(self):
        np.testing.assert_array_equal(rutile.std_types, np.array([14, 14, 8, 8, 8, 8], dtype=np.int32))

    def test_std_positions(self):
        np.testing.assert_array_equal(rutile.std_positions, rutile_dataset['std_positions'])

    def test_std_mapping_to_primitive(self):
        np.testing.assert_array_equal(rutile.std_mapping_to_primitive, np.array([0, 1, 2, 3, 4, 5], dtype=np.int32))

    def test_pointgroup(self):
        self.assertEqual(rutile.pointgroup, '4/mmm')
