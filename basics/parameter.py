#!/usr/bin/env python3
# created at Nov 20, 2017 8:18 PM by Qi Zhang

from collections import namedtuple
from typing import *

Namelist = namedtuple('Namelist', ['name', 'keys'])


def _str_to_QE_str(x: str):
    return x


def _float_to_QE_str(x: float):
    return str(x)


def _bool_to_QE_str(x: bool) -> str:
    if x:
        return '.true.'
    elif not x:
        return '.false.'
    else:
        raise TypeError('Input type unknown!')


def _int_to_QE_str(x: int) -> str:
    return str(x)


class Parameter:
    def __init__(self, name: str, value: Union[str, int, float, bool]):
        """
        Generate an `Parameter` object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param value: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        self._name = name
        self._value = value
        self._type = type(value)

    @property
    def name(self) -> str:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return: the name in INPUTPH_names of the instance
        """
        return self._name

    @property
    def value(self) -> str:
        """
        Raw value given by the user, should be a string.

        :return: the raw value read from input
        """
        return self._value

    @value.setter
    def value(self, new_value: str) -> None:
        self._value = new_value

    @property
    def type(self) -> Type[Union[str, int, float, bool]]:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return:
        """
        return self._type

    def to_QE_str(self) -> str:
        recipe: Dict[Type[Union[str, int, float, bool]], Callable] = {str: _str_to_QE_str,
                                                                      int: _int_to_QE_str,
                                                                      float: _float_to_QE_str,
                                                                      bool: _bool_to_QE_str}
        return recipe[self.type](self.value)

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "{0}: {1}".format(self.name, self.value)

    __repr__ = __str__
