#!/usr/bin/env python3
# created at Nov 20, 2017 8:18 PM by Qi Zhang

from typing import *


class Param:
    def __init__(self, name: str, raw_val: str, namedtuples: Dict[str, NamedTuple]):
        """
        Generate an Param object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param raw_val: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        param: NamedTuple = namedtuples[name]
        self._name = name
        self._raw_value = raw_val
        self._default_value = param[1]  # `default_value` of `param`
        self._type = param[2]  # `type` of `param`

    @property
    def name(self) -> str:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return: the name in INPUTPH_card of the instance
        """
        return self._name

    @property
    def default_value(self):
        """
        The default value given by Quantum ESPRESSO. Read-only property.

        :return: the default value given by Quantum ESPRESSO
        """
        return self._default_value

    @property
    def raw_value(self) -> str:
        """
        Raw value given by the user, should be a string.

        :return: the raw value read from input
        """
        return self._raw_value

    @raw_value.setter
    def raw_value(self, val: str) -> None:
        self._raw_value = val

    @property
    def type(self) -> str:
        """
        The default type defined by Quantum ESPRESSO. Read-only property.

        :return: a string given by Quantum ESPRESSO manual
        """
        return self._type

    @property
    def value(self):
        """
        Convert the raw value given by user to its default type.

        :return: the exact value, with exact type given by Quantum ESPRESSO
        """
        return eval("{0}('{1}')".format(self.type, self.raw_value))

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return '{0}: {1}'.format(self.name, self.value)

    __repr__ = __str__
