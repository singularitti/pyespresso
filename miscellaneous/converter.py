#!/usr/bin/env python3
# created at Aug 7, 2017 2:22 PM by Qi Zhang

from typing import *

import numpy as np


class UnitConverter:
    """
    This is a basic unit converter. It will be handful if you implement further classes based on this.
    """

    @staticmethod
    def simplest_converter(scale: Dict[str, float], num: Union[int, float], from_unit: str, to_unit: str) -> float:
        """
        This is the basic logic of unit-conversion.
        Input a number and 2 different units, it will automatically transform the number from the first to the second
        unit.

        :param scale: a dictionary that has units (string) as one key, and a scaling number from a unit to SI unit as
            the corresponding value.
        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)


class LengthConverter(UnitConverter):
    """
    A converter which converts one-dimensional length.
    """

    def __init__(self):
        """
        Define some constants as scaling numbers.
        """
        self.bohr_radius = 5.2917721067e-11

    def simple_converter(self, num: Union[float, int], from_unit: str, to_unit: str) -> float:
        """
        This function converts the submit first to meter, then converts it in desired unit.

        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        scale: Dict[str, float] = dict.fromkeys(['m', 'meter', 'metre', 'SI'], 1)  # Meter is the standard unit.
        scale.update(dict.fromkeys(['a', 'angstrom', 'A'], 1e-10))
        scale.update(dict.fromkeys(['cm', 'centimeter'], 1e-2))
        scale.update(dict.fromkeys(['nm', 'nanometer'], 1e-9))
        scale.update(dict.fromkeys(['b', 'bohr', 'au', 'atomic'], self.bohr_radius))
        return self.simplest_converter(scale, num, from_unit, to_unit)


class VolumeConverter(UnitConverter):
    """
    A converter which converts three-dimensional volume.
    """

    def __init__(self):
        """
        Define some constants as scaling numbers.
        """
        self.bohr_radius = 5.2917721067e-11

    def simple_converter(self, num, from_unit='a3', to_unit='b3') -> float:
        """
        This function converts the submit first to cubic meter, then converts it in desired unit.

        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        scale: Dict[str, float] = dict.fromkeys(['m3', 'cubicmeter'], 1)
        scale.update(dict.fromkeys(['cm3', 'cubiccentimeter'], 1e-6))
        scale.update(dict.fromkeys(['nm3', 'cubicnanometer'], 1e-27))
        scale.update(dict.fromkeys(['a3', 'cubicangstrom'], 1e-30))
        scale.update(dict.fromkeys(['b3', 'cubicbohr', 'au3', 'cubicatomic'], self.bohr_radius ** 3))
        return self.simplest_converter(scale, num, from_unit, to_unit)


class EnergyConverter(UnitConverter):
    """
    A converter converts energy in different units.
    """

    def __init__(self):
        """
        Define some constants as scaling numbers.
        """
        self.hartree_energy = 4.359744650e-18
        self.electron_volt = 1.602176565e-19
        self.boltzmann_const = 1.38064852e-23
        self.freq_to_joule = 1.98630e-23
        self.hertz_to_joule = 6.62561e-34

    def simple_converter(self, num: Union[float, int], from_unit: str, to_unit: str):
        """
        This function converts the submit first to Joule, then converts it to desired unit.

        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        scale: Dict[str, float] = dict.fromkeys(['J', 'Joule', 'SI'], 1)  # Joule is the standard unit.
        scale.update(dict.fromkeys(['h', 'ha', 'hartree'], self.hartree_energy))
        scale.update(dict.fromkeys(['ev', 'eV', 'electronvolt'], self.electron_volt))
        scale.update(dict.fromkeys(['ry', 'rydberg'], self.hartree_energy / 2))
        scale.update(dict.fromkeys(['K'], self.boltzmann_const))
        scale.update(dict.fromkeys(['cm-1'], self.freq_to_joule))
        scale.update(dict.fromkeys(['hz', 'hertz'], self.hertz_to_joule))
        return self.simplest_converter(scale, num, from_unit, to_unit)


class PressureConverter(UnitConverter):
    """
    A converter converts pressure in different units.
    """

    def simple_converter(self, num: Union[float, int], from_unit: str, to_unit: str) -> float:
        """
        This function converts the submit first to pascal, then converts it in desired unit.

        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        scale: Dict[str, float] = dict.fromkeys(['pa', 'pascal', 'Pa', 'SI'], 1)  # Pa is the standard unit.
        scale.update(dict.fromkeys(['gpa', 'GPa'], 1e9))
        scale.update(dict.fromkeys(['bar'], 1e5))
        scale.update(dict.fromkeys(['mbar', 'Mbar'], 1e11))
        return self.simplest_converter(scale, num, from_unit, to_unit)


class MoleConverter(UnitConverter):
    """
    A converter than convert if there is an unit of form x/mole.
    """

    def __init__(self):
        """
        Define some constants as scaling numbers.
        """
        self.avogadro_const = 6.022140857e23

    def converter(self, num: Union[float, int], from_unit: str, to_unit: str, option: str) -> float:
        """
        Sometimes we meet a quantity with xx/mol, this converts it to the physical quantity of only one such quantity.
        Supported options: {'length', 'volume', 'energy', 'pressure'}

        :param option: str
        :param num: the number to be converted
        :param from_unit: the unit to be converted from
        :param to_unit: the unit to be converted to
        :return: the resulted number
        """
        options = {'l': LengthConverter, 'v': VolumeConverter, 'e': EnergyConverter,
                   'p': PressureConverter}
        return options[option]().simple_converter(num / self.avogadro_const, from_unit, to_unit)


def call_simple_converter(physical_quantity: str, numeric: Union[int, float, List, np.ndarray], from_unit: str,
                          to_unit: str) -> Union[float, List[float], np.ndarray]:
    """
    This is a polymorphism function, it can apply to any object that defined simple_converter method.

    :param physical_quantity: The physical quantity that you want to converter between 2 of its units.
    :param numeric: a float or a list of floats, or a numpy array, in terms of the first unit
    :param from_unit: unit you want to convert from
    :param to_unit: nit you want to convert to
    :return: the resulting value in the second unit
    """
    pool = {'l': LengthConverter, 'v': VolumeConverter, 'e': EnergyConverter,
            'p': PressureConverter}
    converter = pool[physical_quantity]()

    if isinstance(numeric, float) or isinstance(numeric, int):
        return converter.simple_converter(numeric, from_unit, to_unit)
    elif isinstance(numeric, list):
        return converter.simple_converter(np.array(numeric), from_unit, to_unit).tolist()
    elif isinstance(numeric, np.ndarray):
        return converter.simple_converter(numeric, from_unit, to_unit)
    else:
        raise TypeError('Unknown data type' + str(type(numeric)))
