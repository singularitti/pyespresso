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
from typing import Union, List, Tuple, MutableMapping

import addict
from lazy_property import LazyProperty, LazyWritableProperty

from pyque.default_configuerations.namelist import *
from pyque.util.strings import string_to_general_float

# ========================================= What can be exported? =========================================
__all__ = ['is_namelist', 'LazyNamelist', 'NamelistDict', 'builtin_to_qe_string', 'qe_string_to_builtin',
           'DEFAULT_CONTROL_NAMELIST', 'DEFAULT_SYSTEM_NAMELIST', 'DEFAULT_ELECTRONS_NAMELIST', 'DEFAULT_IONS_NAMELIST',
           'DEFAULT_CELL_NAMELIST', 'DEFAULT_INPUTPH_NAMELIST', 'CONTROLNamelistParameter', 'SYSTEMNamelistParameter',
           'ELECTRONSNamelistParameter', 'IONSNamelistParameter', 'CELLNamelistParameter', 'INPUTPHNamelistParameter']

# ========================================= type alias =========================================
NamelistParameterValue = Union[str, int, float, bool]


# ========================================= These are some useful functions. =========================================
def is_namelist(obj: object):
    """
    If an object *obj* is an instance of a class, which is the subclass of abstract base class ``NamelistABC``, then the
    object can be regarded as a namelist.

    :param obj: The object to be examined.
    :return: Whether *obj* is a namelist or not.

    .. doctest::

        >>> is_namelist(DEFAULT_CONTROL_NAMELIST)
        True
    """
    if issubclass(obj.__class__, NamelistABC):
        return True
    return False


def is_namelist_parameter(obj: object):
    """
    If an object *obj* is an instance of a class, which is the subclass of abstract base class ``NamelistParameterABC``,
    then the object can be regarded as a namelist.

    :param obj: The object to be examined.
    :return: Whether *obj* is a namelist parameter or not.

    .. doctest::

        >>> is_namelist_parameter(CONTROLNamelistParameter('calculation', 'scf'))
        True
        >>> is_namelist_parameter(('calculation', 'scf', 'str', 'scf', 'CONTROL'))
        False
    """
    if issubclass(obj.__class__, NamelistParameterABC):
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
    def typed_parameters(self) -> MutableMapping[str, Tuple[Union[str, int, float, bool], str]]:
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
            object.__setattr__(self, name, self[name].__class__(name, value))
        except ValueError:
            raise KeyError("The name '{0}' you set does not already exist!".format(name))

    def beautify(self):
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


# =================================== instances of the crucial class ``Namelist`` ===================================
# These namelists are constants, do not change them.
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

    @property
    @abstractmethod
    def in_namelist(self):
        pass

    @in_namelist.setter
    @abstractmethod
    def in_namelist(self, new_in_namelist):
        pass


class NamelistParameterGeneric(NamelistParameterABC):
    def __init__(self, name: str, value: NamelistParameterValue):
        """
        Generate an `NamelistParameter` object, which stores user given name, and value.

        :param name:
        :param value:
        """
        self.__name = name
        self.__value = value
        self.__default_type = None
        self.__default_value = None
        self.__in_namelist = None

    @property
    def name(self) -> str:
        return self.__name

    @property
    def default_type(self) -> str:
        return self.__default_type

    @default_type.setter
    def default_type(self, new_default_type: str) -> None:
        if new_default_type not in {'int', 'float', 'bool', 'str'}:
            raise ValueError("Unknown default type {0} for name {1} is given!".format(new_default_type, self.__name))
        self.__default_type = new_default_type

    @property
    def value(self) -> NamelistParameterValue:
        return self.__value

    @value.setter
    def value(self, new_value: NamelistParameterValue) -> None:
        self.__value = qe_string_to_builtin(new_value, self.default_type)

    @property
    def default_value(self) -> NamelistParameterValue:
        return self.__default_value

    @default_value.setter
    def default_value(self, new_default_value: NamelistParameterValue) -> None:
        if not isinstance(new_default_value, (int, float, bool, str)):
            raise TypeError("Unknown default value {0} for name {1} is given!".format(new_default_value, self.__name))
        self.__default_value = new_default_value

    @property
    def in_namelist(self):
        return self.__in_namelist

    @in_namelist.setter
    def in_namelist(self, new_in_namelist: str) -> None:
        if new_in_namelist not in {'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'INPUTPH'}:
            raise ValueError("Unknown namelist caption '{0}' is given!".format(new_in_namelist))
        self.__in_namelist = new_in_namelist

    def to_qe_str(self) -> str:
        return builtin_to_qe_string(self.value)

    def to_dict(self):
        return {'name': self.name, 'value': self.value, 'default_type': self.default_type}

    def __str__(self) -> str:
        """
        Show `self.name` and `self.value`.

        :return: a string showing `self.name` and `self.value`
        """
        return "The parameter is '{0}', with value: {1}, and default type: '{2}'.".format(self.name, self.value,
                                                                                          self.default_type)

    def __repr__(self) -> str:
        return ' '.join(map(str, (self.name, self.value, self.default_type, self.default_value, self.in_namelist)))


class NamelistParameter(NamelistParameterGeneric):
    def __init__(self, default_namelist: DefaultNamelist, name: str, value: str):
        super(NamelistParameter, self).__init__(name, value)
        if not is_namelist(default_namelist):
            raise ValueError("The 'default_namelist' given is not a namelist!")
        try:
            default_namelist.names[name]
        except KeyError:
            raise ValueError("The name '{0}' is not a valid name for namelist '{1}'!".format(name, type(self).__name__))
        self.default_value, self.default_type = default_namelist.typed_parameters[name]
        self.in_namelist = default_namelist.caption


class CONTROLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CONTROL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(CONTROLNamelistParameter, self).__init__(DEFAULT_CONTROL_NAMELIST, name, value)


class SYSTEMNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'SYSTEM' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(SYSTEMNamelistParameter, self).__init__(DEFAULT_SYSTEM_NAMELIST, name, value)


class ELECTRONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'ELECTRONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(ELECTRONSNamelistParameter, self).__init__(DEFAULT_ELECTRONS_NAMELIST, name, value)


class IONSNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'IONS' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(IONSNamelistParameter, self).__init__(DEFAULT_IONS_NAMELIST, name, value)


class CELLNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'CELL' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(CELLNamelistParameter, self).__init__(DEFAULT_CELL_NAMELIST, name, value)


class INPUTPHNamelistParameter(NamelistParameter):
    """
    To build a parameter for 'INPUTPH' namelist, you only need to specify the name of the parameter, and its value.
    The type of it will be automatically recognized by its name.
    """

    def __init__(self, name: str, value: str):
        super(INPUTPHNamelistParameter, self).__init__(DEFAULT_INPUTPH_NAMELIST, name, value)


class RawNamelistParameterGeneric(NamelistParameterABC):
    """
    This is used only in parser step.
    """

    def __init__(self, name: str, value: str):
        self.__name = name
        self.__value = value
        self.__in_namelist = None

    def name(self):
        return self.__name

    @property
    def value(self):
        return self.__value

    @value.setter
    def value(self, new_value: str):
        if not isinstance(new_value, str):
            raise TypeError("The type of this value should be 'str'!")
        self.__value = new_value

    @property
    def in_namelist(self):
        return self.__in_namelist

    @in_namelist.setter
    def in_namelist(self, new_in_namelist: str) -> None:
        if new_in_namelist not in {'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'INPUTPH'}:
            raise ValueError("Unknown namelist caption '{0}' is given!".format(new_in_namelist))
        self.__in_namelist = new_in_namelist

    @abstractmethod
    def eval(self):
        pass


class RawNamelistParameter(RawNamelistParameterGeneric):
    def __init__(self, default_namelist: DefaultNamelist, name: str, value: str):
        super(RawNamelistParameter, self).__init__(name, value)
        if not is_namelist(default_namelist):
            raise ValueError("The 'default_namelist' given is not a namelist!")
        try:
            default_namelist.names[name]
        except KeyError:
            raise ValueError("The name '{0}' is not a valid name for namelist '{1}'!".format(name, self.__class__))
        self.in_namelist = default_namelist.caption

    @abstractmethod
    def eval(self):
        pass


class RawCONTROLNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super(RawCONTROLNamelistParameter, self).__init__(DEFAULT_CONTROL_NAMELIST, name, value)

    def eval(self):
        return CONTROLNamelistParameter(self.__name, self.__value)


class RawSYSTEMNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super(RawSYSTEMNamelistParameter, self).__init__(DEFAULT_SYSTEM_NAMELIST, name, value)

    def eval(self):
        return SYSTEMNamelistParameter(self.__name, self.__value)


class RawELECTRONSNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super(RawELECTRONSNamelistParameter, self).__init__(DEFAULT_ELECTRONS_NAMELIST, name, value)

    def eval(self):
        return ELECTRONSNamelistParameter(self.__name, self.__value)


class RawIONSNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super().__init__(DEFAULT_IONS_NAMELIST, name, value)

    def eval(self):
        return IONSNamelistParameter(self.__name, self.__value)


class RawCELLNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super().__init__(DEFAULT_CELL_NAMELIST, name, value)

    def eval(self):
        return CELLNamelistParameter(self.__name, self.__value)


class RawINPUTPHNamelistParameter(RawNamelistParameter):
    def __init__(self, name: str, value: str):
        super().__init__(DEFAULT_INPUTPH_NAMELIST, name, value)

    def eval(self):
        return INPUTPHNamelistParameter(self.__name, self.__value)

# ================================================= end of this block =================================================
