#!/usr/bin/env python3
"""
:mod:`card` -- 
========================================

.. module card
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from lazy_property import LazyWritableProperty


class LazyCard(LazyWritableProperty):
    pass
    # def __set__(self, instance, value):
    #     """
    #
    #     :param instance: Here is a PWscfInput instance.
    #     :param value: should be a dict or NamelistDict instance.
    #     :return:
    #     """
    #     print(value)
    #     print(CardDict({self.cache_name: value}))
    #     super().__set__(instance, CardDict({self.cache_name: value}))
