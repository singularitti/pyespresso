#!/usr/bin/env python3
"""
:mod:`namelist` -- title
========================================

.. module namelist
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import re
from abc import ABC, abstractmethod
from typing import Union, Dict

import addict
from lazy_property import LazyWritableProperty

from pyque.default_configuerations.namelist import *
from pyque.tools.strings import string_to_general_float, all_string_like

# ========================================= What can be exported? =========================================
__all__ = ['LazyNamelist', 'NamelistDict', 'builtin_to_qe_string', 'qe_string_to_builtin', 'NamelistVariable',
           'CONTROLNamelistVariable', 'SYSTEMNamelistVariable',
           'ELECTRONSNamelistVariable', 'IONSNamelistVariable', 'CELLNamelistVariable', 'INPUTPHNamelistVariable']

# ========================================= type alias =========================================
NamelistVariableValue = Union[str, int, float, bool]


# ========================================= These are some useful functions. =========================================
def is_namelist_variable(obj: object):
    """
    If an object *obj* is an instance of a class, which is the subclass of abstract base class ``NamelistVariableABC``,
    then the object can be regarded as a namelist.

    :param obj: The object to be examined.
    :return: Whether *obj* is a namelist parameter or not.

    .. doctest::

        >>> is_namelist_variable(CONTROLNamelistVariable('calculation', 'scf'))
        True
    """
    if issubclass(obj.__class__, NamelistVariableABC):
        return True
    return False


def builtin_to_qe_string(obj: object):
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
        raise TypeError("Invalid type '{0}' is given!".format(type(obj).__name__))


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
        raise TypeError("Invalid type '{0}' is given!".format(type(obj).__name__))


# ========================================= define some crucial class =========================================
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
            object.__setattr__(self, name, self[name].__class__(name, value))
        except ValueError:
            raise KeyError("The name '{0}' you set does not already exist!".format(name))

    @LazyWritableProperty
    def _group_name(self):
        pass

    def eval(self):
        d = dict()
        for k, v in self.items():
            d.update({k: v.__class__(v.name, v.value)})
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


# ========================================= define other crucial classes =========================================
class NamelistVariableABC(ABC):
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

    @property
    @abstractmethod
    def in_namelist(self):
        pass

    @in_namelist.setter
    @abstractmethod
    def in_namelist(self, new_in_namelist):
        pass


class NamelistVariable(NamelistVariableABC):
    def __init__(self, group_name: str, name: str, value: NamelistVariableValue):
        # Type checking
        if not all_string_like({group_name, name}):
            raise TypeError("Arguments *in_namelist* and *name* must both be strings!")
        if not isinstance(value, (int, float, bool, str)):
            raise TypeError("Argument *value* is of wrong type '{0}'!".format(type(value)))
        # Validity checking
        if group_name not in {'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'INPUTPH'}:
            raise ValueError("Unknown namelist caption '{0}' is given!".format(group_name))
        # Some names have numbers as their labels, like 'celldm(1)'. So we need to separate them.
        if '(' in name:
            # Only take the part before the first '('
            name_prefix = re.match("(\w+)\(?(\d*)\)?", name, flags=re.IGNORECASE).group(1)
        else:
            name_prefix = name
        typed_namelist = TYPED_NAMELISTS[group_name]
        if name_prefix in typed_namelist:
            self.__name = name
            self.__default_type = typed_namelist[name_prefix]
            self.__in_namelist = group_name
        else:
            raise ValueError("The name '{0}' is not valid in namelist '{1}'!".format(name, group_name))
        self.__value = value

    @property
    def name(self) -> str:
        return self.__name

    @property
    def default_type(self) -> str:
        return self.__default_type

    @property
    def value(self) -> NamelistVariableValue:
        return self.__value

    @value.setter
    def value(self, new_value: NamelistVariableValue) -> None:
        self.__value = qe_string_to_builtin(new_value, self.default_type)

    @property
    def in_namelist(self) -> str:
        return self.__in_namelist

    def to_qe_str(self) -> str:
        return builtin_to_qe_string(self.value)

    def to_dict(self) -> Dict[str, Union[str, NamelistVariableValue]]:
        return {'name': self.name, 'value': self.value, 'default_type': self.default_type}

    def __eq__(self, other: object) -> bool:
        attributes = ['name', 'default_type', 'value', 'in_namelist']
        return all(getattr(self, attr) == getattr(other, attr) for attr in attributes)

    def __ne__(self, other) -> bool:
        attributes = ['name', 'default_type', 'value', 'in_namelist']
        return any(getattr(self, attr) != getattr(other, attr) for attr in attributes)


class CONTROLNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(CONTROLNamelistVariable, self).__init__('CONTROL', name, value)


class SYSTEMNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(SYSTEMNamelistVariable, self).__init__('SYSTEM', name, value)


class ELECTRONSNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(ELECTRONSNamelistVariable, self).__init__('ELECTRONS', name, value)


class IONSNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(IONSNamelistVariable, self).__init__('IONS', name, value)


class CELLNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(CELLNamelistVariable, self).__init__('CELL', name, value)


class INPUTPHNamelistVariable(NamelistVariable):
    """
    To build a parameter for 'INPUTPH' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(INPUTPHNamelistVariable, self).__init__('INPUTPH', name, value)

# ================================================= end of this block =================================================
