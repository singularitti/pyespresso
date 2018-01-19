#!/usr/bin/env python3
"""
:mod:`cell` -- An adaptive for ``spglib`` library
=================================================

.. module:: atom_mass
   :platform: Unix, Windows, Mac, Linux
   :synopsis: This module integrates spglib for Python API, and built 2 classes: ``SimpleCell`` for simply storing some basic
   information about a cell, and ``Cell`` for much more complex functional extension.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import os
from numbers import Number
from typing import *

import numpy as np
import spglib
from json_tricks import dump
from lazy_property import LazyProperty

# ========================================= What can be exported? =========================================
__all__ = ['is_simple_cell', 'is_cell', 'get_symmetry_from_database', 'get_spacegroup_type',
           'get_hall_number_from_symmetry', 'Cell', 'SimpleCell', 'simple_cell_to_cell', 'print_cell']

# ================================= These are some type aliases or type definitions. =================================
CellInitialValue = TypeVar('CellInitialValue', list, np.ndarray, float, bool, int, None)

# =================================== These are some exports of spglib functions. ===================================
get_symmetry_from_database: Callable[[int], Dict[str, np.ndarray]] = spglib.get_symmetry_from_database

get_spacegroup_type: Callable[[int], Dict[str, Union[str, int]]] = spglib.get_spacegroup_type

get_hall_number_from_symmetry: Callable[[np.ndarray, np.ndarray, float], int] = spglib.get_hall_number_from_symmetry


# ========================================= These are some useful functions. =========================================
def is_simple_cell(obj: object) -> bool:
    """
    If an object is an instance of ``SimpleCell``, then it can be regarded as a simple cell.

    :param obj: The object to be checked.
    :return: Whether it is a simple cell.
    """
    if isinstance(obj, SimpleCell):
        return True
    else:
        return False


def is_cell(obj: object) -> bool:
    """
    If an object is an instance of ``Cell``, then it can be regarded as a cell.

    :param obj: The object to be checked.
    :return: Whether it is a cell.
    """
    if isinstance(obj, Cell):
        return True
    else:
        return False


def _cell_initial_values_equal(a: CellInitialValue, b: CellInitialValue) -> bool:
    """
    For the initial values for 2 cells, i.e., 'lattice', 'positions', 'numbers', 'magmoms', '_symprec',
    '_angle_tolerance', '_symbol_type', '_eps', and '_silent', compare them according to their different data
    structures.

    :param a: An initial value for one cell.
    :param b: The corresponding initial value for another cell.
    :return: Whether a and b are equal.
    """
    if any(type(x) == np.ndarray for x in (a, b)):
        return np.array_equal(a, b)
    else:
        return a == b


# ============================== Code block for some validating functions ==============================
def _validate_lattice(lattice: np.ndarray):
    if not lattice.shape == (3, 3):
        raise ValueError('The shape of ``lattice`` is not 3x3!')
    if not issubclass(lattice.dtype.type, (Number, np.number)):
        raise TypeError('The lattice parameters are not made of numbers!')


def _validate_positions(positions: np.ndarray, atoms_num: int):
    shape = positions.shape
    if not shape == (atoms_num, 3):
        raise ValueError('The shape of ``positions`` is not Nx3!')
    if not issubclass(positions.dtype.type, (Number, np.number)):
        raise TypeError('The positions are not made of numbers!')


def _validate_numbers(numbers: np.ndarray, atoms_num: int):
    shape = numbers.shape
    if not shape == (atoms_num,):
        raise ValueError('The shape of ``numbers`` is not of shape Nx1!')
    if not issubclass(numbers.dtype.type, (Number, np.number)):
        raise TypeError('The ``numbers`` is not made of numbers!')


def _validate_magmoms(magmoms: np.ndarray, atoms_num: int):
    shape = magmoms.shape
    if not shape == (atoms_num,):
        raise ValueError('The shape of ``magmoms`` is not of shape Nx1!')
    if not issubclass(magmoms.dtype.type, (Number, np.number)):
        raise TypeError('The ``magmoms`` is not made of numbers!')


# ====================================== The followings are classes definitions. ======================================
class Cell:
    def __init__(self, lattice: Union[list, np.ndarray], positions: Union[list, np.ndarray],
                 numbers: Union[list, np.ndarray], *args):
        """
        More detailed documentation see `here <https://atztogo.github.io/spglib/python-spglib.html#python-spglib/>`_.

        :param lattice: Lattice parameters are given by a $3x3$ matrix with floating point values, where a,b,c
            are given as rows, which results in the transpose of the definition for C-API.
        :param positions: Fractional atomic positions positions are given by a $Nx3$ matrix with floating point values,
            where $N$ is the number of atoms.
        :param numbers: Numbers to distinguish atomic species. Numbers are given by a list of N integers.
        :param args: Only 1 optional positional argument is allowed, which is ``magmoms`` in ``spglib`` library.
            The collinear polarizations magmoms only work with ``get_symmetry`` and are given as a list of $N$ floating
            point values.
        """

        # Why we need to convert all ``lattice``, ``positions``, ``numbers`` and ``magmoms`` (if given) to numpy arrays?
        # Because in ``refine_cell``, ``find_primitive`` and ``standardize_cell`` we generate new cells, and they are
        # already numpy arrays according to spglib. So if we want to this will make people confused if the resulted
        # attributes change type. And will cause ambiguity when comparing cells (using ``__eq__`` and ``__ne__``).

        self.lattice: np.ndarray = np.array(lattice)
        self.positions: np.ndarray = np.array(positions)
        self.numbers: np.ndarray = np.array(numbers)
        self.atoms_num = len(self.numbers)

        _validate_lattice(self.lattice)
        _validate_positions(self.positions, self.atoms_num)
        _validate_numbers(self.numbers, self.atoms_num)

        if len(args) == 0:
            self.magmoms = None
            self.cell: Tuple[np.ndarray, ...] = (self.lattice, self.positions, self.numbers)
        elif len(args) == 1:
            self.magmoms: np.ndarray = np.array(args[0])
            _validate_magmoms(self.magmoms, self.atoms_num)
            self.cell: Tuple[np.ndarray, ...] = (self.lattice, self.positions, self.numbers, self.magmoms)
        else:
            raise TypeError(
                'Only 1 optional positional argument ``magmoms`` is allowed, but {0} are given!'.format(len(args)))

        self._symprec: float = 1e-5
        self._angle_tolerance: float = -1.0
        self._symbol_type: int = 0
        self._eps: float = 1e-5
        self._silent = False

    # ============================== Code block for some "hidden" values ==============================
    @property
    def symprec(self):
        return self._symprec

    @symprec.setter
    def symprec(self, new_symprec: float):
        if not self.silent:
            print('Be careful, setting symprec may affect the values of some attributes!')
        self._symprec = new_symprec

    @property
    def angle_tolerance(self):
        return self._angle_tolerance

    @angle_tolerance.setter
    def angle_tolerance(self, new_angle_tolerance: float):
        if not self.silent:
            print('Be careful, setting angle_tolerance may affect the values of some attributes!')
        self._angle_tolerance = new_angle_tolerance

    @property
    def symbol_type(self):
        return self._symbol_type

    @symbol_type.setter
    def symbol_type(self, new_symbol_type: int):
        if type(new_symbol_type) is not int:
            raise TypeError('Symbol type should be an integer, but {0} given!'.format(type(new_symbol_type)))
        if not self.silent:
            print('Be careful, setting symbol_type may affect the values of some attributes!')
        self._symbol_type = new_symbol_type

    @property
    def eps(self):
        return self._eps

    @eps.setter
    def eps(self, new_eps: float):
        if not self.silent:
            print('Be careful, setting eps may affect the values of some attributes!')
        self._eps = new_eps

    @property
    def silent(self):
        return self._silent

    @silent.setter
    def silent(self, new_silent: bool):
        self._silent = new_silent

    # ========================================= end =========================================

    def get_spacegroup(self, symprec: float = 1e-5, angle_tolerance: float = -1.0, symbol_type: int = 0) -> str:
        """

        :param symprec:
        :param angle_tolerance:
        :param symbol_type:
        :return:
        """
        return spglib.get_spacegroup(self.cell, symprec, angle_tolerance, symbol_type)

    @LazyProperty
    def spacegroup(self) -> str:
        """
        A property as a shorthand of ``get_spacegroup``.

        :return:
        """
        return self.get_spacegroup(symprec=self.symprec, angle_tolerance=self.angle_tolerance,
                                   symbol_type=self.symbol_type)

    def get_symmetry(self, symprec: float = 1e-5, angle_tolerance=-1.0) -> Dict[str, np.array]:
        return spglib.get_symmetry(self.cell, symprec, angle_tolerance)

    @LazyProperty
    def symmetry(self) -> Dict[str, np.array]:
        return self.get_symmetry(symprec=self.symprec, angle_tolerance=self.angle_tolerance)

    # ============================== Code block for structure optimization ==============================
    def refine_cell(self, symprec: float = 1e-5, angle_tolerance: float = -1.0) -> Optional['Cell']:
        """
        This is just a wrapper for ``spglib.refine_cell``.

        :param symprec:
        :param angle_tolerance:
        :return:
        """
        search_result = spglib.refine_cell(self.cell, symprec, angle_tolerance)
        if not search_result:  # If ``search_result`` is ``None``
            print('Refine cell failed!')
        else:
            lattice, scaled_positions, numbers = search_result
            return Cell(lattice, scaled_positions, numbers)

    def find_primitive(self, symprec: float = 1e-5, angle_tolerance: float = -1.0) -> Optional['Cell']:
        """
        This is just a wrapper for ``spglib.find_primitive``.

        :param symprec:
        :param angle_tolerance:
        :return:
        """
        search_result = spglib.find_primitive(self.cell, symprec, angle_tolerance)
        if not search_result:  # If ``search_result`` is ``None``
            print('Find primitive cell failed!')
        else:
            lattice, scaled_positions, numbers = search_result
            return Cell(lattice, scaled_positions, numbers)

    def standardize_cell(self, to_primitive: bool = False, no_idealize: bool = False, symprec: float = 1e-5,
                         angle_tolerance=-1.0) -> Optional['Cell']:
        """
        This is just a wrapper for ``spglib.standardize_cell``.

        :param to_primitive:
        :param no_idealize:
        :param symprec:
        :param angle_tolerance:
        :return:
        """
        search_result = spglib.standardize_cell(self.cell, to_primitive, no_idealize, symprec, angle_tolerance)
        if not search_result:  # If ``search_result`` is ``None``
            print('Standardized crystal structure search failed!')
        else:
            lattice, scaled_positions, numbers = search_result
            return Cell(lattice, scaled_positions, numbers)

    # ========================================= end =========================================

    # ============================== Code block for symmetry dataset ==============================
    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1.0, hall_number=0) -> Optional[dict]:
        """

        :param symprec:
        :param angle_tolerance:
        :param hall_number:
        :return:
        """
        d: dict = spglib.get_symmetry_dataset(self.cell, symprec, angle_tolerance, hall_number)
        if not d:  # If ``ds`` is ``None``
            print('The search failed!')
        else:
            return d

    @LazyProperty
    def symmetry_dataset(self) -> Optional[Dict[str, Union[int, np.ndarray, str, List[str]]]]:
        """
        A property as a shorthand of ``get_symmetry_dataset`` method.
        Note only default arguments are used! They are: ``symprec=1e-5, angle_tolerance=-1.0, hall_number=0``.
        If you want a more flexible control of parameters, you can either:

        * Use ``get_symmetry_dataset`` method directly and give arguments you want.
        * change ``self.symprec``, ``self.angle_tolerance``, ``self.hall_number``. But be careful, if you change
          these values, the returned symmetry dataset may be different from previous one, and thus cause some
          other attributes different.

        :return: Symmetry dataset.
        """
        return self.get_symmetry_dataset(symprec=self.symprec, angle_tolerance=self.angle_tolerance)

    @LazyProperty
    def number(self) -> int:
        """
        This is just a wrapper for ``symmetry_dataset['number']``.

        :return:
        """
        return self.symmetry_dataset['number']

    @LazyProperty
    def hall_number(self) -> int:
        """
        This is just a wrapper for ``symmetry_dataset['hall_number']``.

        :return:
        """
        return self.symmetry_dataset['hall_number']

    @LazyProperty
    def international(self) -> str:
        """
        This is just a wrapper for ``symmetry_dataset['international']``.

        :return:
        """
        return self.symmetry_dataset['international']

    international_short = international

    @LazyProperty
    def hall(self) -> str:
        """
        This is just a wrapper for ``symmetry_dataset['hall']``.

        :return:
        """
        return self.symmetry_dataset['hall']

    hall_symbol = hall

    @LazyProperty
    def choice(self) -> str:
        """
        This is just a wrapper for ``symmetry_dataset['choice']``.

        :return:
        """
        return self.symmetry_dataset['choice']

    @LazyProperty
    def transformation_matrix(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['transformation_matrix']``.

        :return:
        """
        return self.symmetry_dataset['transformation_matrix']

    @LazyProperty
    def origin_shift(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['origin_shift']``.

        :return:
        """
        return self.symmetry_dataset['origin_shift']

    @LazyProperty
    def rotations(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['rotations']``.

        :return:
        """
        return self.symmetry_dataset['rotations']

    @LazyProperty
    def translations(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['translations']``.

        :return:
        """
        return self.symmetry_dataset['translations']

    @LazyProperty
    def wyckoffs(self) -> List[str]:
        """
        This is just a wrapper for ``symmetry_dataset['wyckoffs']``.

        :return:
        """
        return self.symmetry_dataset['wyckoffs']

    @LazyProperty
    def equivalent_atoms(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['equivalent_atoms']``.

        :return:
        """
        return self.symmetry_dataset['equivalent_atoms']

    @LazyProperty
    def mapping_to_primitive(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['mapping_to_primitive']``.

        :return:
        """
        return self.symmetry_dataset['mapping_to_primitive']

    @LazyProperty
    def std_lattice(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['std_lattice']``.

        :return:
        """
        return self.symmetry_dataset['std_lattice']

    @LazyProperty
    def std_types(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['std_types']``.

        :return:
        """
        return self.symmetry_dataset['std_types']

    @LazyProperty
    def std_positions(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['std_positions']``.

        :return:
        """
        return self.symmetry_dataset['std_positions']

    @LazyProperty
    def std_mapping_to_primitive(self) -> np.ndarray:
        """
        This is just a wrapper for ``symmetry_dataset['std_mapping_to_primitive']``.

        :return:
        """
        return self.symmetry_dataset['std_mapping_to_primitive']

    @LazyProperty
    def pointgroup(self) -> str:
        """
        This is just a wrapper for ``symmetry_dataset['pointgroup']``.

        :return:
        """
        return self.symmetry_dataset['pointgroup']

    # ========================================= end =========================================

    # ============================== Code block for reduction methods ==============================
    def niggli_reduce(self, eps: float = 1e-5):
        """
        This is just a wrapper for ``spglib.niggli_lattice``.

        :param eps:
        :return:
        """
        search_result = spglib.niggli_reduce(self.lattice, eps)
        if not search_result:  # If ``search_result`` is ``None``
            print('Niggli reduction search failed!')
        else:
            niggli_lattice, = search_result
            return niggli_lattice

    @LazyProperty
    def niggli_lattice(self):
        """
        This is just a property as a shorthand for ``self.niggli_reduce`` method.

        :return:
        """
        return self.niggli_reduce(eps=self.eps)

    def delaunay_reduce(self, eps=1e-5):
        """
        This is just a wrapper for ``spglib.delaunay_reduce``.

        :param eps:
        :return:
        """
        search_result = spglib.delaunay_reduce(self.lattice, eps)
        if not search_result:  # If ``search_result`` is ``None``
            print('Delaunay reduction search failed!')
        else:
            delaunay_lattice, = search_result
            return delaunay_lattice

    @LazyProperty
    def delaunay_lattice(self):
        """
        This is just a property as a shorthand for ``self.delaunay_reduce`` method.

        :return:
        """
        return self.delaunay_reduce(eps=self.eps)

    def get_ir_reciprocal_mesh(self, mesh: List[int], is_shift: List[int] = [0, 0, 0]):
        """
        This is just a wrapper for ``spglib.get_ir_reciprocal_mesh``.

        :param mesh:
        :param is_shift:
        :return:
        """
        search_result = spglib.get_ir_reciprocal_mesh(mesh, self.cell, is_shift)
        if not search_result:  # If ``search_result`` is ``None``
            print('Delaunay reduction search failed!')
        else:
            mapping, grid = search_result
            return mapping, grid

    # ========================================= end =========================================

    # ============================== Code block for spacegroup type ==============================
    @LazyProperty
    def spacegroup_type(self) -> Dict[str, Union[str, int]]:
        """

        :return:
        """
        # Note that ``self.hall_number`` comes from ``self.symmetry_dataset['hall_number']``, so if
        # ``self.symmetry_dataset`` is changed by accident, unwanted error may arise!
        return get_spacegroup_type(self.hall_number)

    @LazyProperty
    def international_full(self) -> str:
        """
        This is just a wrapper for ``spacegroup_type['international_full']``.

        :return:
        """
        return self.spacegroup_type['international_full']

    @LazyProperty
    def schoenflies(self) -> str:
        """
        This is just a wrapper for ``spacegroup_type['schoenflies']``.

        :return:
        """
        return self.spacegroup_type['schoenflies']

    @LazyProperty
    def pointgroup_schoenflies(self) -> str:
        """
        This is just a wrapper for ``spacegroup_type['pointgroup_schoenflies']``.

        :return:
        """
        return self.spacegroup_type['pointgroup_schoenflies']

    @LazyProperty
    def pointgroup_international(self) -> str:
        """
        This is just a wrapper for ``spacegroup_type['pointgroup_international']``.

        :return:
        """
        return self.spacegroup_type['pointgroup_international']

    @LazyProperty
    def arithmetic_crystal_class_number(self) -> int:
        """
        This is just a wrapper for ``spacegroup_type['arithmetic_crystal_class_number']``.

        :return:
        """
        return self.spacegroup_type['arithmetic_crystal_class_number']

    @LazyProperty
    def arithmetic_crystal_class_symbol(self) -> str:
        """
        This is just a wrapper for ``spacegroup_type['arithmetic_crystal_class_symbol']``.

        :return:
        """
        return self.spacegroup_type['arithmetic_crystal_class_symbol']

    # ========================================= end =========================================

    # ============================== Code block for special methods ==============================
    def __str__(self) -> str:
        return \
            """
            lattice: {0}
            positions: {1}
            numbers: {2}
            """.format(self.lattice.tolist(), self.positions.tolist(), self.numbers.tolist())

    def __repr__(self):
        return {'lattice': self.lattice, 'positions': self.positions, 'numbers': self.numbers}

    def __eq__(self, other: object) -> bool:
        """
        Two cells are equal if all of their ``lattice``, ``positions``, ``numbers``, ``magmoms``, ``symprec``,
        ``angle_tolerance``, ``symbol_type``, ``eps`` and ``silent`` attributes are equal.
        Why only compare these attributes? Because they uniquely define the cell.

        :param other: Should be another cell, or else raise an error.
        :return: Whether two cells are equal.
        """
        if is_cell(other):
            attrs = ['lattice', 'positions', 'numbers', 'magmoms', '_symprec', '_angle_tolerance', '_symbol_type',
                     '_eps', '_silent']
            return all(_cell_initial_values_equal(self.__dict__[x], other.__dict__[x]) for x in attrs)
        else:
            raise TypeError('{0} is not a cell type!'.format(other))

    def __ne__(self, other: object) -> bool:
        """
        Two cells are not equal if any of their ``lattice``, ``positions``, ``numbers``, ``magmoms``, ``symprec``,
        ``angle_tolerance``, ``symbol_type``, ``eps`` and ``silent`` attributes is not equal.
        Why only compare these attributes? Because they uniquely define the cell.

        :param other: Should be another cell, or else raise an error.
        :return: Whether two cells are not equal.
        """
        if is_cell(other):
            attrs = ['lattice', 'positions', 'numbers', 'magmoms', '_symprec', '_angle_tolerance', '_symbol_type',
                     '_eps', '_silent']
            for x in attrs:  # Find the first attribute which are not equal for the 2 cells, then exit
                if not _cell_initial_values_equal(self.__dict__[x], other.__dict__[x]):
                    print('Attribute {0} for the 2 cells are not equal!'.format(x))
                    return True
            else:
                return False
        else:
            raise TypeError('{0} is not a cell type!'.format(other))

    def to_json(self, output_name: str, output_path: str = ''):
        if not output_name.endswith('.json'):
            output_name: str = output_name + '.json'
        if not output_path:
            file_path: str = os.path.join(os.getcwd(), output_name)
        else:
            file_path: str = os.path.join(output_path, output_name)

        with open(file_path, 'w') as f:
            dump(self.__repr__(), f)

    # ========================================= end =========================================


class SimpleCell:
    def __init__(self, lattice: Union[list, np.ndarray], positions: Union[list, np.ndarray],
                 numbers: Union[list, np.ndarray]):
        self.lattice: np.ndarray = np.array(lattice)
        self.positions: np.ndarray = np.array(positions)
        self.numbers: np.ndarray = np.array(numbers)
        self.atoms_num = len(self.numbers)

        _validate_lattice(self.lattice)
        _validate_positions(self.positions, self.atoms_num)
        _validate_numbers(self.numbers, self.atoms_num)

    def to_cell(self) -> Cell:
        return Cell(self.lattice, self.positions, self.numbers)

    def __str__(self) -> str:
        return \
            """
            lattice: {0}
            positions: {1}
            numbers: {2}
            """.format(self.lattice.tolist(), self.positions.tolist(), self.numbers.tolist())

    def __repr__(self):
        return {'lattice': self.lattice, 'positions': self.positions, 'numbers': self.numbers}

    def __eq__(self, other: object) -> bool:
        """
        Two simple cells are equal if all of their ``lattice``, ``positions``, ``numbers`` attributes are equal.

        :param other: Should be another simple cell, or else raise an error.
        :return: Whether two cells are equal.
        """
        if is_simple_cell(other):
            attrs = ['lattice', 'positions', 'numbers']
            return all(np.array_equal(self.__dict__[x], other.__dict__[x]) for x in attrs)
        else:
            raise TypeError('{0} is not a simple cell type!'.format(other))

    def __ne__(self, other: object) -> bool:
        """
        Two simple cells are not equal if any of their ``lattice``, ``positions``, ``numbers`` attributes is not equal.

        :param other: Should be another simple cell, or else raise an error.
        :return: Whether two cells are not equal.
        """
        if is_simple_cell(other):
            attrs = ['lattice', 'positions', 'numbers']
            for x in attrs:  # Find the first attribute which are not equal for the 2 cells, then exit
                if not np.array_equal(self.__dict__[x], other.__dict__[x]):
                    print('Attribute {0} for the 2 simple cells are not equal!'.format(x))
                    return True
            else:
                return False
        else:
            raise TypeError('{0} is not a cell type!'.format(other))

    def to_json(self, output_name: str, output_path: str = ''):
        if not output_name.endswith('.json'):
            output_name: str = output_name + '.json'
        if not output_path:
            file_path: str = os.path.join(os.getcwd(), output_name)
        else:
            file_path: str = os.path.join(output_path, output_name)

        with open(file_path, 'w') as f:
            dump(self.__repr__(), f)


def print_cell(cell: Union[Cell, SimpleCell]) -> None:
    """
    Pretty print a (simple) cell.

    :param cell: The (simple) cell to be printed.
    :return: None.
    """
    if is_cell(cell) or is_simple_cell(cell):
        try:
            from beeprint import pp
            pp(cell)
        except ModuleNotFoundError:
            print(cell)
    else:
        raise TypeError('{0} is not a cell!'.format(cell))


def simple_cell_to_cell(sc: SimpleCell) -> Optional[Cell]:
    """
    Suppose you have a simple cell, and you want to convert it to a more complex cell, then use this method.

    :param sc: The simple cell to be converted.
    :return: The cell which is converted result.
    """
    if is_simple_cell(sc):
        return sc.to_cell()
    else:
        raise TypeError('{0} is not a simple cell!'.format(sc))
