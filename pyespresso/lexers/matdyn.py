#!/usr/bin/env python3
"""
:mod:`mod` -- 
========================================

.. module mod
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from lazy_property import LazyWritableProperty


class FreqLexer:
    @LazyWritableProperty
    def text_stream(self):
        pass

