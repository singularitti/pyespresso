#!/usr/bin/env python3
# created at Aug 7, 2017 2:22 PM by Nil-Zil


class UnitConverter:
    """
    This is a unit converter.
    Input a number and 2 different units,
    it will automatically transform the first in terms of the second one.
    """

    def __init__(self):
        """
        Define some SI unit constants.
        """
        self.bohr_radius = 5.2917721067e-11
        self.hartree_energy = 4.359744650e-18
        self.electron_volt = 1.602176565e-19
        self.boltzmann_const = 1.38064852e-23
        self.avogadro_const = 6.022140857e23

    def length_converter(self, num, from_unit='a', to_unit='b'):
        """
        This function converts the input first to meter, then converts it to desired unit.
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        scale = dict.fromkeys(['m', 'meter', 'metre', 'SI'],
                              1)  # Meter is the standard unit.
        scale.update(dict.fromkeys(['a', 'angstrom', 'A'], 1e-10))
        scale.update(dict.fromkeys(['cm', 'centimeter'], 1e-2))
        scale.update(dict.fromkeys(['nm', 'nanometer'], 1e-9))
        scale.update(dict.fromkeys(
            ['b', 'bohr', 'au', 'atomic'], self.bohr_radius))
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)

    def volume_converter(self, num, from_unit='a3', to_unit='b3'):
        """
        This function converts the input first to cubic meter, then converts it to desired unit.
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        scale = dict.fromkeys(['m3', 'cubicmeter'], 1)
        scale.update(dict.fromkeys(['cm3', 'cubiccentimeter'], 1e-6))
        scale.update(dict.fromkeys(['nm3', 'cubicnanometer'], 1e-27))
        scale.update(dict.fromkeys(['a3', 'cubicangstrom'], 1e-30))
        scale.update(dict.fromkeys(
            ['b3', 'cubicbohr', 'au3', 'cubicatomic'], self.bohr_radius ** 3))
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)

    def energy_converter(self, num, from_unit='ha', to_unit='ry'):
        """
        This function converts the input first to Joule, then converts it to desired unit.
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        scale = dict.fromkeys(['J', 'Joule', 'SI'],
                              1)  # Joule is the standard unit.
        scale.update(dict.fromkeys(
            ['h', 'ha', 'hartree'], self.hartree_energy))
        scale.update(dict.fromkeys(
            ['ev', 'eV', 'electronvolt'], self.electron_volt))
        scale.update(dict.fromkeys(['ry', 'rydberg'], self.hartree_energy / 2))
        scale.update(dict.fromkeys(['K'], self.boltzmann_const))
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)

    @staticmethod
    def pressure_converter(num, from_unit='gpa', to_unit='mbar'):
        """
        This function converts the input first to pascal, then converts it to desired unit.
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        scale = dict.fromkeys(['pa', 'pascal', 'Pa', 'SI'],
                              1)  # Pa is the standard unit.
        scale.update(dict.fromkeys(['gpa', 'GPa'], 1e9))
        scale.update(dict.fromkeys(['bar'], 1e5))
        scale.update(dict.fromkeys(['mbar', 'Mbar'], 1e11))
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)

    def mole_converter(self, option, num, from_unit, to_unit):
        """
        Sometimes we meet a quantity with xx/mol, this converts it to the physical quantity of only one such quantity.
        Supported options: {'length', 'volume', 'energy', 'pressure'}
        :param option: str
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        options = {'length': self.length_converter, 'volume': self.volume_converter, 'energy': self.energy_converter,
                   'pressure': self.pressure_converter}
        return options[option](num / self.avogadro_const, from_unit, to_unit)


if __name__ == '__main__':
    uc = UnitConverter()
    # print(uc.length_converter(7.4268, 'b', 'a'))
    # print(uc.volume_converter(11.176, 'a3', 'b3'))
    print(uc.energy_converter(2000, 'K', 'ry'))
    # print(uc.pressure_converter(3, 'mbar', 'gpa'))
    # print(uc.mole_converter('volume', 6.73, 'cm3', 'au3'))
