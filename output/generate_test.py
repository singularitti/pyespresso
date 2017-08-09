#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 23, 2017 10:26 PM by Nil-Zil
"""
Quickly generate testing cases for plotting and checking.
"""

import re


class GenerateTest:
    def __init__(self, prefix: str, infix: list, suffix=".dat"):
        """
        :param prefix:
        :param infix:
        :param suffix:
        """
        self.prefix = prefix
        self.infix = infix
        self.suffix = suffix

    def generate_filelist(self) -> list:
        """

        :return: list(str)
        """
        return [self.prefix + inf + self.suffix for inf in self.infix]

    def generate_legend(self) -> list:
        """
        To use this function, keep each of your varaiable separated by '_' or '-'.
        :return: list
        """
        return [' '.join(re.split(r'[-_]', leg)) for leg in self.infix]
