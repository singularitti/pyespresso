#!/usr/bin/env python3
# created at Nov 20, 2017 8:18 PM by Qi Zhang

from typing import *


class Param:
    def __init__(self, name: str, raw_value: Union[str, int, float, bool]):
        """
        Generate an Param object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param raw_value: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        self._name = name
        self._raw_value = raw_value

    @property
    def name(self) -> str:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return: the name in INPUTPH_names of the instance
        """
        return self._name

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
    def value(self):
        """
        Convert the raw value given by user to its default type.

        :return: the exact value, with exact type given by Quantum ESPRESSO
        """
        return self.type(self.raw_value)

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "{0}: {1}".format(self.name, self.value)

    __repr__ = __str__
