#!/usr/bin/env python3
# created at Nov 20, 2017 8:18 PM by Qi Zhang

from collections import defaultdict
from typing import DefaultDict, Type, Union, Tuple, List, Callable, Dict

# Type alias
DefaultParameters = DefaultDict[str, Tuple[Union[str, int, float, bool], Type[Union[str, int, float, bool]]]]


def _str_to_QE_str(x: str) -> str:
    return x


def _float_to_QE_str(x: float) -> str:
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
    def __init__(self, name: str, value: Union[str, int, float, bool], value_type: Type[Union[str, int, float, bool]]):
        """
        Generate an `Parameter` object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param value: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        self._name = name
        self._value = value
        self._type = value_type
        self._default_value = None
        self._in_namelist = None

    @property
    def name(self) -> str:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return: the name in INPUTPH_names of the instance
        """
        return self._name

    @property
    def value(self) -> Union[str, int, float, bool]:
        """
        Raw value given by the user, should be a string.

        :return: the raw value read from input
        """
        return self._type(self._value)

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

    @property
    def default_value(self) -> Union[str, int, float, bool]:
        """
        The default value given by Quantum ESPRESSO.
        Read-only property.

        :return: the default value given by Quantum ESPRESSO
        """
        return self._default_value

    @property
    def in_namelist(self) -> str:
        """
        Return the namelist the parameter is in.
        Read-only property.

        :return: namelist's name
        """
        return self._in_namelist

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
        return "The parameter is '{0}', value: {1}, with type: {2}".format(self.name, self.value, self.type.__name__)

    def __repr__(self):
        return "'{0}': {1}".format(self.name, self.value)


class Namelist:
    def __init__(self, namelist_name: str, names: List[str], default_values: List):
        self.__name__ = namelist_name
        self._names = names
        self._default_values = default_values

    @property
    def names(self) -> List[str]:
        """
        This is a list of names for a namelist.

        :return:
        """
        return self._names

    @property
    def default_values(self) -> List[Union[str, int, float, bool]]:
        """
        This is a list of QE default values for a namelist.

        :return:
        """
        return self._default_values

    @property
    def value_types(self) -> List[Type[Union[str, int, float, bool]]]:
        """
        This is a list of types of each default value for a namelist.

        :return:
        """
        return [type(x) for x in self._default_values]

    @property
    def default_parameters(self) -> DefaultParameters:
        """
        This is a `DefaultDict` of `(QE default value, types of each default value)` tuples for a namelist, with those
        names to be its keys.

        :return:
        """
        default_parameters: DefaultParameters = defaultdict(tuple)
        default_values: List = self._default_values
        value_types: List = self.value_types
        for i, k in enumerate(self._names):
            default_parameters[k] = (default_values[i], value_types[i])
        return default_parameters
