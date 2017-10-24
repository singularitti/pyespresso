#!/usr/bin/env python3
# created at Jul 23, 2017 10:26 PM by Qi Zhang
"""
Quickly generate testing cases for plot.
"""

from read_file.elasticity import *


class GenerateTest:
    def __init__(self):
        self.sr = SimpleReader()

    def generate_from_ls(self, test_filename) -> list:
        return self.sr.read_each_line(test_filename)

    @staticmethod
    def generate_from_composition(prefix: str, infix: list, suffix=""):
        return [prefix + inf + suffix for inf in infix]
