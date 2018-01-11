#!/usr/bin/env python3
"""
:mod:`parameter` -- title
========================================

.. module parameter
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from collections import namedtuple
from typing import Union, Type, Dict, Callable, NamedTuple

from pyque.meta.namelist import Namelist, CONTROL_NAMELIST, SYSTEM_NAMELIST, ELECTRONS_NAMELIST, IONS_NAMELIST, \
    CELL_NAMELIST, INPUTPH_NAMELIST

# ========================================= What can be exported? =========================================
__all__ = ['CONTROLParameter', 'SYSTEMParameter', 'ELECTRONSParameter', 'IONSParameter', 'CELLParameter',
           'INPUTPHParameter']

# ================================= These are some type aliases or type definitions. =================================
SimpleParameter = NamedTuple('SimpleParameter', [('name', str), ('value', Union[str, int, float, bool]), ('type', str)])

# ========================================= define useful data structures =========================================
SimpleParameter: SimpleParameter = namedtuple('SimpleParameter', ['name', 'value', 'type'])


def _str_to_qe_str(x: str) -> str:
    return x


def _float_to_qe_str(x: float) -> str:
    return str(x)


def _bool_to_qe_str(x: bool) -> str:
    if x:
        return '.true.'
    elif not x:
        return '.false.'
    else:
        raise TypeError('Input type unknown!')


def _int_to_qe_str(x: int) -> str:
    return str(x)


class ParameterGeneric:
    def __init__(self, name: str, value: Union[str, int, float, bool], value_type: Type[Union[str, int, float, bool]]):
        """
        Generate an `ParameterGeneric` object, which stores user given name, and value.

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
    def type(self) -> str:
        """
        Read-only property, cannot be changed once the instance is generated.

        :return:
        """
        return self._type.__name__

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

    def to_qe_str(self) -> str:
        recipe: Dict[str, Callable] = {'str': _str_to_qe_str,
                                       'int': _int_to_qe_str,
                                       'float': _float_to_qe_str,
                                       'bool': _bool_to_qe_str}
        return recipe[self.type.__name__](self.value)

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "The parameter is '{0}', with value: {1}, and type: {2}.".format(self._name, self.value,
                                                                                self.type)

    def __repr__(self):
        return str(SimpleParameter(self._name, self.value, self.type))


class _Parameter(ParameterGeneric):
    """
    This is a prototype for a parameter.
    You only need to provide the name of your parameter, the value of it, and the parameter ``dict`` it belongs to.
    """

    def __init__(self, name: str, value: str, namelist: Namelist):
        self._default_value, value_type = namelist.default_parameters[name]
        super().__init__(name, value, value_type)
        self._in_namelist = namelist.__name__


class CONTROLParameter(_Parameter):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CONTROL_NAMELIST)


class SYSTEMParameter(_Parameter):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, SYSTEM_NAMELIST)


class ELECTRONSParameter(_Parameter):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, ELECTRONS_NAMELIST)


class IONSParameter(_Parameter):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, IONS_NAMELIST)


class CELLParameter(_Parameter):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CELL_NAMELIST)


class INPUTPHParameter(_Parameter):
    def __init__(self, name, value):
        super().__init__(name, value, INPUTPH_NAMELIST)
