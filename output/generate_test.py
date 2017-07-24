#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 23, 2017 10:26 PM by Nil-Zil
"""
Quickly generate testing cases for plotting and checking.
"""


class GenerateTest(object):
    def __init__(self, prefix: str, infix: list, suffix=".dat"):
        """
        :param prefix:
        :param infix:
        :param suffix:
        """
        self.prefix = prefix
        self.infix = infix
        self.suffix = suffix

    def generate_filelist(self):
        """

        :return: list(str)
        """
        return [self.prefix + infix + self.suffix for infix in self.infix]


