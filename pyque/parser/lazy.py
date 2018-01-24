#!/usr/bin/env python3
"""
:mod:`lazy` -- 
========================================

.. module lazy
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import re
from typing import *

from lazy_property import LazyWritableProperty


class NamelistParser(LazyWritableProperty):
    @staticmethod
    def get(instance):
        return instance

    def for_namelist(self) -> str:
        return re.match("parse_(\w+)_namelist", self.cache_name).group(1).upper()

    def data(self) -> Dict[str, str]:
        return self.method()

