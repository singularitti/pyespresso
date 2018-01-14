#!/usr/bin/env python3
"""
:mod:`namelist` -- title
========================================

.. module namelist
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Union, List, DefaultDict, Tuple

import addict
from lazy_property import LazyProperty, LazyWritableProperty

from pyque.default_configuerations.namelist import *
from pyque.util.strings import string_to_general_float

# ========================================= What can be exported? =========================================
__all__ = ['is_namelist', 'LazyNamelist', 'NamelistDict', 'builtin_to_qe_string', 'qe_string_to_builtin']


# ========================================= These are some useful functions. =========================================
def is_namelist(obj: object):
    if issubclass(type(obj), NamelistABC):
        return True
    return False


def is_namelist_parameter(obj: object):
    if issubclass(type(obj), NamelistParameterABC):
        return True
    return False


def builtin_to_qe_string(obj: Union[int, float, bool, str]):
    """
    An instance of a builtin type can be converted to a Quantum ESPRESSO legal string.

    :param obj: Legal type of *obj* can be one of ``int``, ``str``, ``bool``, ``float``, cannot be anything else. If you
        have, for example, a list of builtins, map this function on them.
    :return: A Quantum ESPRESSO string.
    """
    if isinstance(obj, bool):
        # You must check this at first, since ``bool`` is a subclass of ``int``, if you check this below
        # the ``elif`` clause, then ``isinstance(True, bool)`` will be ``True``.
        if obj:  # If *obj* is ``True``.
            return '.true.'
        return '.false.'  # If *obj* is ``False``.
    elif isinstance(obj, (int, float, str)):
        return str(obj)
    else:
        raise TypeError("Illegal type '{0}' is given!".format(type(obj).__name__))


def qe_string_to_builtin(obj: Union[int, float, bool, str], desired_type: str):
    """
    A Quantum ESPRESSO legal string can be converted to an instance of a Python builtin type.

    :param obj: Legal type of *obj* can be one of ``int``, ``str``, ``bool``, ``float``, cannot be anything else. This
        maybe a little bit confusing because the name of this function is *qe_string_to_builtin*. However, it is for
        compatibility consideration because you may already have a an instance of a Python builtin type.
    :param desired_type: Should be a string indicating the type you want the *obj* to be converted to.
    :return: An instance of a Python builtin type, converted from *obj*.
    """
    if desired_type not in {'int', 'str', 'bool', 'float'}:
        raise ValueError("Unknown *desired_type* '{0}' is given!".format(desired_type))
    if isinstance(obj, (int, float, bool)):  # In case that your data are already builtins.
        if type(obj).__name__ == desired_type:
            return obj
        raise TypeError(
            "The type of '{0}' ({1}) is inconsistent with *desired_type* {2}!".format(obj, type(obj).__name__,
                                                                                      desired_type))
    elif isinstance(obj, str):
        if desired_type == 'bool':
            try:
                return {'.true.': True, '.false.': False}[obj]
            except KeyError:
                raise ValueError("The string '{0}' cannot be converted to a bool!".format(obj))
        elif desired_type == 'int':
            try:
                return int(obj)
            except ValueError:
                raise ValueError('The string {0} cannot be converted to an integer!'.format(obj))
        elif desired_type == 'float':
            try:
                return string_to_general_float(obj)
            except ValueError:
                raise ValueError('The string {0} cannot be converted to a float!'.format(obj))
        elif desired_type == 'str':
            return obj
    else:  # Here *obj* is not ``int``, ``str``, ``bool``, ``float``.
        raise TypeError("Illegal type '{0}' is given!".format(type(obj).__name__))


# ========================================= define some crucial class =========================================
class NamelistABC(ABC):
    @abstractmethod
    def names(self):
        pass

    @abstractmethod
    def values(self):
        pass

    @abstractmethod
    def value_types(self):
        pass


class DefaultNamelist(NamelistABC):
    """
    This will build a constant instance which is a constant. So the names of this class instances should be all
    capitalized.
    """

    def __init__(self, caption: str, odict: OrderedDict):
        self.__caption = caption
        self.parameters = odict

    @LazyProperty
    def caption(self):
        return self.__caption

    @LazyProperty
    def names(self) -> List[str]:
        """
        This is a list of names for a namelist. Read-only attribute.

        :return:
        """
        return list(self.parameters.keys())

    @LazyProperty
    def values(self) -> List[Union[str, int, float, bool]]:
        """
        This is a list of QE default values for a namelist. Read-only attribute.

        :return:
        """
        return list(self.parameters.values())

    @LazyProperty
    def value_types(self) -> List[str]:
        """
        This is a list of types of each default value for a namelist. Read-only attribute.

        :return:
        """
        return [type(_).__name__ for _ in self.values]

    @LazyProperty
    def typed_parameters(self) -> DefaultDict[str, Tuple[Union[str, int, float, bool], str]]:
        """
        This is a ``DefaultDict`` of ``(QE default value, types of each default value)`` tuples for a namelist,
        with those names to be its keys. Read-only attribute.

        :return:
        """
        return OrderedDict(zip(self.names, zip(self.values, self.value_types)))


class LazyNamelist(LazyWritableProperty):
    def __set__(self, instance, value):
        """

        :param instance: For here it can be a PWscfInput instance.
        :param value: should be a dict or NamelistDict instance.
        :return:
        """
        if isinstance(value, dict):
            super().__set__(instance, NamelistDict(value))
        else:
            raise ValueError('You should set it to be a dict!')


class NamelistDict(addict.Dict):
    def __setattr__(self, name, value):
        try:
            object.__setattr__(self, name, eval(type(self[name]).__name__)(name, value))
        except NameError:
            raise NameError("The name '{0}' you set does not already exist!".format(name))

    def beautify(self):
        d = dict()
        for k, v in self.items():
            d.update({k: type(v)(v.name, v.value)})
        return NamelistDict(d)

    def to_dict(self):
        d = dict()
        for k, v in self.items():
            d.update({k: v.to_dict()})
        return d

    def to_text(self):
        lines = []
        for k, v in self.items():
            lines.append("{0} = {1}".format(k, builtin_to_qe_string(v.value)))
        return '\n'.join(lines)


# ================================================= end of this block =================================================


# =================================== instances of the crucial class ``Namelist`` ===================================
DEFAULT_CONTROL_NAMELIST: DefaultNamelist = DefaultNamelist('CONTROL', CONTROL_NAMELIST_DICT)

DEFAULT_SYSTEM_NAMELIST: DefaultNamelist = DefaultNamelist('SYSTEM', SYSTEM_NAMELIST_DICT)

DEFAULT_ELECTRONS_NAMELIST: DefaultNamelist = DefaultNamelist('ELECTRONS', ELECTRONS_NAMELIST_DICT)

DEFAULT_IONS_NAMELIST: DefaultNamelist = DefaultNamelist('IONS', IONS_NAMELIST_DICT)

DEFAULT_CELL_NAMELIST: DefaultNamelist = DefaultNamelist('CELL', CELL_NAMELIST_DICT)

DEFAULT_INPUTPH_NAMELIST = DefaultNamelist('INPUTPH', INPUTPH_NAMELIST_DICT)


# ================================================= end of this block =================================================


# ========================================= define other crucial classes =========================================
class NamelistParameterABC(ABC):
    @abstractmethod
    def name(self):
        pass

    @property
    @abstractmethod
    def value(self):
        pass

    @value.setter
    @abstractmethod
    def value(self, new_value):
        pass

    @abstractmethod
    def default_type(self):
        pass

    @abstractmethod
    def in_namelist(self):
        pass


class NamelistParameter(NamelistParameterABC):
    def __init__(self, name: str, value: Union[str, int, float, bool]):
        """
        Generate an `NamelistParameterMeta` object, which stores user given name, and value.

        :param name: the name given by user in the card
        :param value: a raw value given by user, has to be a string, it will be converted into exact type defined by
            `self.type` parameter.
        """
        self.__name = name
        self.__value = value

    @LazyProperty
    def name(self) -> str:
        return self.__name

    @LazyWritableProperty
    def _default_type(self) -> str:
        pass

    @LazyProperty
    def default_type(self) -> str:
        return self._default_type

    @property
    def value(self) -> Union[str, int, float, bool]:
        return self.__value

    @value.setter
    def value(self, new_value: Union[str, int, float, bool]) -> None:
        self.__value = qe_string_to_builtin(new_value, self.default_type)

    @LazyWritableProperty
    def _default_value(self) -> Union[str, int, float, bool]:
        pass

    @LazyProperty
    def default_value(self) -> Union[str, int, float, bool]:
        return self._default_value

    @LazyWritableProperty
    def _in_namelist(self) -> str:
        pass

    @LazyProperty
    def in_namelist(self) -> str:
        return self._in_namelist

    def to_qe_str(self) -> str:
        return builtin_to_qe_string(self.value)

    def to_dict(self):
        return {'name': self.name, 'real_value': self.value, 'intrinsic_type': self.default_type}

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "The parameter is '{0}', with value: {1}, and default type: '{2}'.".format(self.name, self.value,
                                                                                          self.default_type)

    def __repr__(self) -> str:
        return str((self.name, self.value, self.default_type, self.default_value, self.in_namelist))


class CONTROLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_CONTROL_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_CONTROL_NAMELIST.caption


class SYSTEMNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_SYSTEM_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_SYSTEM_NAMELIST.caption


class ELECTRONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_ELECTRONS_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_ELECTRONS_NAMELIST.caption


class IONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_IONS_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_IONS_NAMELIST.caption


class CELLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_CELL_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_CELL_NAMELIST.caption


class INPUTPHNamelistParameter(NamelistParameter):
    def __init__(self, name, value):
        super().__init__(name, value)
        self._default_value, self._default_type = DEFAULT_INPUTPH_NAMELIST.typed_parameters[name]
        self._in_namelist = DEFAULT_INPUTPH_NAMELIST.caption
# ================================================= end of this block =================================================
