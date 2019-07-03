#!/usr/bin/env python3
"""
:mod:`console` -- 
========================================

.. module console
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import code
from collections import OrderedDict
from typing import *
from lazy_property import LazyProperty

from pyespresso.core.namelist import DEFAULT_CONTROL_NAMELIST, CONTROLNamelistVariable


class PWscfInteractiveConsole:
    @LazyProperty
    def banner(self):
        return """Welcome to '{0}'! Now input things you want. 
        And we will automatically generate an input file for you!""".format(self.__class__)

    @LazyProperty
    def exit_message(self):
        return """Exiting... Thank you for using '{0}'! Have a nice one!""".format(self.__class__)

    # def raw_input(self, prompt=""):
    #     self.write(prompt)
    #     line = self.instream.readline()
    #     if line:
    #         return line.rstrip("\n")
    #     raise EOFError()

    def generate_control_namelist(self):
        cs = code.InteractiveConsole()
        odict: MutableMapping[str, Tuple[Union[str, int, float, bool], str]] = DEFAULT_CONTROL_NAMELIST.typed_variables
        user_dict = OrderedDict()
        for name, (value, value_type) in odict.items():
            raw_val = cs.raw_input(prompt="Please input the value you want for name '{0}', \
                the type of it is '{1}', the default value of it is {2}:".format(name, value_type, value))
            user_dict[name] = CONTROLNamelistVariable(name, raw_val)
        return user_dict

    # def generate_system_namelist(self):
    # return self.raw_input()
