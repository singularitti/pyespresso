#!/usr/bin/env python3
# created at Nov 28, 2017 12:52 AM by Qi Zhang

from typing import *

import numpy as np
import spglib

from basics.lazy import CachedProperty

get_symmetry_from_database: Callable = spglib.get_symmetry_from_database

get_spacegroup_type: Callable = spglib.get_spacegroup_type

get_hall_number_from_symmetry: Callable = spglib.get_hall_number_from_symmetry


class Cell:
    def __init__(self, lattice, positions, numbers, *args):
        """
        More detailed documentation see [here](https://atztogo.github.io/spglib/python-spglib.html#python-spglib).

        :param lattice: Lattice parameters lattice are given by a $3x3$ matrix with floating point values, where a,b,c
            are given as rows, which results in the transpose of the definition for C-API.
        :param positions: Fractional atomic positions positions are given by a $Nx3$ matrix with floating point values,
            where $N$ is the number of atoms.
        :param numbers: Numbers to distinguish atomic species numbers are given by a list of N integers.
        :param args: Only 1 optional positional argument is allowed, which is `magmoms` in `spglib` library.
            The collinear polarizations magmoms only work with `get_symmetry` and are given as a list of $N$ floating
            point values.
        """
        self.lattice = lattice
        self.positions = positions
        self.numbers = numbers
        if len(args) == 0:
            self.cell = (self.lattice, self.positions, self.numbers)
        elif len(args) == 1:
            self.magmoms, = args
            self.cell = (self.lattice, self.positions, self.numbers, self.magmoms)
        else:
            raise TypeError(
                'Only 1 optional positional argument `magmoms` is allowed, but {0} given!'.format(len(args)))

    def get_spacegroup(self, symprec: float = 1e-5) -> str:
        return spglib.get_spacegroup(self.cell, symprec)

    def get_symmetry(self, symprec: float = 1e-5) -> Dict[str, np.array]:
        return spglib.get_symmetry(self.cell, symprec)

    def refine_cell(self, symprec: float = 1e-5):
        pass

    def find_primitive(self, symprec: float = 1e-5):
        return self.standardize_cell(to_primitive=True, no_idealize=False, symprec=symprec)

    def standardize_cell(self, to_primitive: bool = False, no_idealize: bool = False, symprec: float = 1e-5) -> \
            Optional[object]:
        """
        Implemented using `spglib.standardize_cell`.

        :param to_primitive:
        :param no_idealize:
        :param symprec:
        :return:
        """
        standardized = spglib.standardize_cell(self.cell, to_primitive, no_idealize, symprec)
        if not standardized:  # If `standardized` is `None`
            print('Standardized crystal structure search failed!')
        else:
            lattice, scaled_positions, numbers = standardized
            return Cell(lattice, scaled_positions, numbers, None)

    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1.0, hall_number=0) -> dict:
        ds = spglib.get_symmetry_dataset(self.cell, symprec, angle_tolerance, hall_number)
        if not ds:  # If `ds` is `None`
            print('The search failed!')
        else:
            return ds

    def get_symmetry_from_database(self):
        hall_number = self.symmetry_database['hall_number']
        return get_symmetry_from_database(hall_number)

    def get_spacegroup_type(self):
        hall_number = self.symmetry_database['hall_number']
        return get_spacegroup_type(hall_number)

    def get_hall_number_from_symmetry(self, symprec=1e-5):
        rotations, translations = self.symmetry_database['rotations'], self.symmetry_database['translations']
        return get_hall_number_from_symmetry(rotations, translations, symprec)

    def niggli_reduce(self, eps: float = 1e-5):
        search_result = spglib.niggli_reduce(self.lattice, eps)
        if not search_result:  # If `search_result` is `None`
            print('Niggli reduction search failed!')
        else:
            niggli_lattice, = search_result
            return niggli_lattice

    def delaunay_reduce(self, eps=1e-5):
        search_result = spglib.delaunay_reduce(self.lattice, eps)
        if not search_result:  # If `search_result` is `None`
            print('Delaunay reduction search failed!')
        else:
            delaunay_lattice, = search_result
            return delaunay_lattice

    def get_ir_reciprocal_mesh(self, mesh: List[int], is_shift: List[int] = [0, 0, 0]):
        search_result = spglib.get_ir_reciprocal_mesh(mesh, self.cell, is_shift)
        if not search_result:  # If `search_result` is `None`
            print('Delaunay reduction search failed!')
        else:
            mapping, grid = search_result
            return mapping, grid

    @CachedProperty
    def niggli_lattice(self):
        return self.niggli_reduce()

    @CachedProperty
    def delaunay_lattice(self):
        return self.delaunay_reduce()

    @CachedProperty
    def spacegroup(self):
        return self.get_spacegroup()

    @CachedProperty
    def symmetry(self):
        return self.get_symmetry()

    @CachedProperty
    def dataset(self):
        return self.get_symmetry_dataset()

    @CachedProperty
    def number(self):
        return self.dataset['number']

    @CachedProperty
    def hall_number(self):
        return self.dataset['hall_number']

    @CachedProperty
    def international(self):
        return self.dataset['international']

    @CachedProperty
    def hall(self):
        return self.dataset['hall']

    @CachedProperty
    def choice(self):
        return self.dataset['choice']

    @CachedProperty
    def transformation_matrix(self):
        return self.dataset['transformation_matrix']

    @CachedProperty
    def origin_shift(self):
        return self.dataset['origin_shift']

    @CachedProperty
    def rotations(self):
        return self.dataset['rotations']

    @CachedProperty
    def translations(self):
        return self.dataset['translations']

    @CachedProperty
    def wyckoffs(self):
        return self.dataset['wyckoffs']

    @CachedProperty
    def equivalent_atoms(self):
        return self.dataset['equivalent_atoms']

    @CachedProperty
    def mapping_to_primitive(self):
        return self.dataset['mapping_to_primitive']

    @CachedProperty
    def std_lattice(self):
        return self.dataset['std_lattice']

    @CachedProperty
    def std_types(self):
        return self.dataset['std_types']

    @CachedProperty
    def std_positions(self):
        return self.dataset['std_positions']

    @CachedProperty
    def std_mapping_to_primitive(self):
        return self.dataset['std_mapping_to_primitive']

    @CachedProperty
    def pointgroup(self):
        return self.dataset['pointgroup']

    def __str__(self):
        return str(self.cell)

    __repr__ = __str__

    def __getattr__(self, item):
        if item == 'symmetry_database':
            self.__dict__['symmetry_database'] = self.get_symmetry_from_database()
            return self.symmetry_database
        else:
            raise AttributeError('Cell does not have attribute {0}!'.format(item))
