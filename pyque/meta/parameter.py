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
from typing import Union, Type, NamedTuple

from lazy_property import LazyProperty, LazyWritableProperty

from pyque.meta.namelist import Namelist, CONTROL_NAMELIST, SYSTEM_NAMELIST, ELECTRONS_NAMELIST, IONS_NAMELIST, \
    CELL_NAMELIST, INPUTPH_NAMELIST
from pyque.lexer.simple import ValueWithComment

# ========================================= What can be exported? =========================================
__all__ = ['CONTROLNamelistParameter', 'SYSTEMNamelistParameter', 'ELECTRONSNamelistParameter', 'IONSNamelistParameter',
           'CELLNamelistParameter', 'INPUTPHNamelistParameter', 'to_qe_str']

# ================================= These are some type aliases or type definitions. =================================
SimpleParameter = NamedTuple('SimpleParameter', [('name', str), ('value', Union[str, int, float, bool]), ('type', str)])

# ========================================= define useful data structures =========================================
SimpleParameter: SimpleParameter = namedtuple('SimpleParameter', ['name', 'value', 'type'])


def to_qe_str(obj):
    if isinstance(obj, (int, float, str)):
        return str(obj)
    elif isinstance(obj, bool):
        if obj:  # If *obj* is ``True``.
            return '.true.'
        return '.false.'  # If *obj* is ``False``.
    # elif isinstance(obj, NamelistParameterMeta):
    #     return to_qe_str(obj.value)
    elif isinstance(obj, ValueWithComment):
        if obj.comment:
            return '   !'.join([to_qe_str(obj.value), obj.comment])
        else:
            return to_qe_str(obj.value)
    raise TypeError('Type {0} given is not legal here!'.format(type(obj).__name__))


class NamelistParameterMeta:
    def __init__(self, name: str, value: Union[str, int, float, bool], value_type: Type[Union[str, int, float, bool]]):
        """
        Generate an `NamelistParameterMeta` object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param value: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        self._name = name
        self._value = value
        self._type = value_type
        self._default_value = None
        self._in_namelist = None

    @LazyProperty
    def name(self) -> str:
        return self._name

    @LazyWritableProperty
    def value(self) -> Union[str, int, float, bool]:
        return self._type(self._value)

    @LazyProperty
    def type(self) -> str:
        return self._type.__name__

    @LazyProperty
    def default_value(self) -> Union[str, int, float, bool]:
        return self._default_value

    @LazyProperty
    def in_namelist(self) -> str:
        return self._in_namelist

    def to_qe_str(self) -> str:
        return to_qe_str(self.value)

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "The parameter is '{0}', with value: {1}, and type: {2}.".format(self._name, self.value, self.type)

    def __repr__(self):
        return str(SimpleParameter(self._name, self.value, self.type))


class NamelistParameter(NamelistParameterMeta):
    """
    This is a prototype for a parameter.
    You only need to provide the name of your parameter, the value of it, and the parameter ``dict`` it belongs to.
    """

    def __init__(self, name: str, value: str, namelist: Namelist):
        self._default_value, value_type = namelist.default_parameters[name]
        super().__init__(name, value, value_type)
        self._in_namelist = namelist.__name__


class CONTROLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CONTROL_NAMELIST)


class SYSTEMNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, SYSTEM_NAMELIST)


class ELECTRONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, ELECTRONS_NAMELIST)


class IONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, IONS_NAMELIST)


class CELLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value, CELL_NAMELIST)


class INPUTPHNamelistParameter(NamelistParameter):
    def __init__(self, name, value):
        super().__init__(name, value, INPUTPH_NAMELIST)
