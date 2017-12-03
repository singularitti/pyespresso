#!/usr/bin/env python3
# created at Dec 2, 2017 2:58 AM by Qi Zhang

import unittest

from basics.cell import *
from basics.unittests.cell_test_data import *

rutile = Cell(*rutile_spglib_cell)
rutile_distorted = Cell(*rutile_distorted_spglib_cell)
simple_rutile = SimpleCell(*rutile_spglib_cell)


class CellTester(unittest.TestCase):
    def test_symprec(self):
        self.assertEqual(rutile.symprec, 1e-5)
        rutile.symprec = 1e-6
        self.assertEqual(rutile.symprec, 1e-6)

    def test_angle_tolerance(self):
        self.assertEqual(rutile.angle_tolerance, -1.0)

    def test_get_spacegroup(self):
        self.assertEqual(rutile.get_spacegroup(), 'P4_2/mnm (136)')

    def test_get_symmetry(self):
        rs = rutile.get_symmetry()
        for k in rutile_symmetry.keys():
            np.testing.assert_array_equal(rs[k], rutile_symmetry[k])

    def test_refine_cell(self):
        refined_cell = rutile_distorted.refine_cell(symprec=1e-1)
        np.testing.assert_array_equal(refined_cell.lattice, np.array([[4., 0., 0.],
                                                                      [0., 4., 0.],
                                                                      [0., 0., 3.]]))
        np.testing.assert_almost_equal(refined_cell.positions, np.array([[0.00000, 0.00000, 0.00000],
                                                                         [0.50000, 0.50000, 0.50000],
                                                                         [0.29995, 0.29995, 0.00000],
                                                                         [0.70005, 0.70005, 0.00000],
                                                                         [0.20005, 0.79995, 0.50000],
                                                                         [0.79995, 0.20005, 0.50000]]))
        np.testing.assert_array_equal(refined_cell.numbers, np.array([14, 14, 8, 8, 8, 8]))

    def test_standardize_cell(self):
        standardized_cell = rutile_distorted.standardize_cell(to_primitive=False, no_idealize=True, symprec=1e-1)
        np.testing.assert_array_equal(standardized_cell.lattice, np.array([[3.97, 0., 0.],
                                                                           [0., 4.03, 0.],
                                                                           [0., 0., 3.]]))
        np.testing.assert_array_equal(standardized_cell.positions, np.array([[0., 0., 0.],
                                                                             [0.5001, 0.5, 0.5],
                                                                             [0.3, 0.3, 0.],
                                                                             [0.7, 0.7, 0.002],
                                                                             [0.2, 0.8, 0.5],
                                                                             [0.8, 0.2, 0.5]]))
        np.testing.assert_array_equal(standardized_cell.numbers, np.array([14, 14, 8, 8, 8, 8]))

    def test_find_primitive(self):
        primitive = rutile_distorted.find_primitive(symprec=1e-1)
        np.testing.assert_array_equal(primitive.lattice, np.array([[4., 0., 0.],
                                                                   [0., 4., 0.],
                                                                   [0., 0., 3.]]))
        np.testing.assert_almost_equal(primitive.positions, np.array([[0.00000, 0.00000, 0.00000],
                                                                      [0.50000, 0.50000, 0.50000],
                                                                      [0.29995, 0.29995, 0.00000],
                                                                      [0.70005, 0.70005, 0.00000],
                                                                      [0.20005, 0.79995, 0.50000],
                                                                      [0.79995, 0.20005, 0.50000]]))
        np.testing.assert_array_equal(primitive.numbers, np.array([14, 14, 8, 8, 8, 8]))

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

    def test_international_short(self):
        self.assertEqual(rutile.international_short, 'P4_2/mnm')

    def test_hall(self):
        self.assertEqual(rutile.hall, '-P 4n 2n')

    def test_hall_symbol(self):
        self.assertEqual(rutile.hall_symbol, '-P 4n 2n')

    def test_choice(self):
        self.assertEqual(rutile.choice, '')

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

    def test_international_full(self):
        self.assertEqual(rutile.international_full, 'P 4_2/m 2_1/n 2/m')

    def test_schoenflies(self):
        self.assertEqual(rutile.schoenflies, 'D4h^14')

    def test_pointgroup_schoenflies(self):
        self.assertEqual(rutile.pointgroup_schoenflies, '4/mmm')

    def test_pointgroup_international(self):
        self.assertEqual(rutile.pointgroup_international, 'D4h')

    def test_arithmetic_crystal_class_number(self):
        print(rutile.arithmetic_crystal_class_number)
        self.assertEqual(rutile.arithmetic_crystal_class_number, 36)

    def test_arithmetic_crystal_class_symbol(self):
        self.assertEqual(rutile.arithmetic_crystal_class_symbol, '4/mmmP')

    @staticmethod
    def test__str__():
        print(rutile)

    def test__eq__(self):
        print(rutile.symprec)
        self.assertFalse(rutile == rutile_distorted)
        self.assertTrue(rutile == rutile)

    def test__ne__(self):
        self.assertTrue(rutile != rutile_distorted)
        self.assertFalse(rutile != rutile)

    def test_is_cell(self):
        self.assertTrue(is_cell(rutile))
        self.assertFalse(any(is_cell(x) for x in [1, 'test', {'lattice': 1, 'positions': 2, 'numbers': 3}]))

    @staticmethod
    def test_print_cell():
        print_cell(rutile)


class SimpleCellTester(unittest.TestCase):
    def test_to_cell(self):
        self.assertTrue(simple_rutile.to_cell() == rutile)
        self.assertTrue(simple_cell_to_cell(simple_rutile) == rutile)

    @staticmethod
    def test__str__():
        print(simple_rutile)

    def test__eq__(self):
        self.assertTrue(simple_rutile == simple_rutile)

    def test__ne__(self):
        self.assertFalse(simple_rutile != simple_rutile)
