#!/usr/bin/env python3
# created at Aug 7, 2017 2:22 PM by Nil-Zil


class UnitConverter(object):
    """
    This is a unit converter.
    """

    def __init__(self):
        """
        Define some SI unit constants.
        """
        self.bohr_radius = 5.2917721067e-11
        self.hartree_energy = 4.359744650e-18
        self.electron_volt = 1.602176565e-19
        self.boltzmann_const = 1.38064852e-23

    def length_converter(self, num, from_unit='a', to_unit='b'):
        """
        This function converts the input first to meter, then converts it to desired unit.
        :param num: float
        :param from_unit: str
        :param to_unit: str
        :return: float
        """
        scale = dict.fromkeys(['m', 'meter', 'metre', 'SI'], 1)  # Meter is the standard unit.
        scale.update(dict.fromkeys(['a', 'angstrom', 'A'], 1e-10))
        scale.update(dict.fromkeys(['nm', 'nanometer'], 1e-9))
        scale.update(dict.fromkeys(['b', 'bohr', 'atomic'], self.bohr_radius))
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
        scale = dict.fromkeys(['J', 'Joule'], 1)  # Joule is the standard unit.
        scale.update(dict.fromkeys(['h', 'ha', 'hartree'], self.hartree_energy))
        scale.update(dict.fromkeys(['ev', 'eV', 'electronvolt'], self.electron_volt))
        scale.update(dict.fromkeys(['ry', 'rydberg'], self.hartree_energy / 2))
        scale.update(dict.fromkeys(['K'], self.boltzmann_const))
        normalized_num = num * scale.get(from_unit)
        return normalized_num / scale.get(to_unit)


if __name__ == '__main__':
    uc = UnitConverter()
    print(uc.length_converter(0.3808, 'nm', 'b'))
    print(uc.energy_converter(1000, 'ev', 'J'))
