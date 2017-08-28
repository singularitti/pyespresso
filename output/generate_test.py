#!/usr/bin/env python3
# created at Jul 23, 2017 10:26 PM by Nil-Zil
"""
Quickly generate testing cases for plotting.
"""

from output.read_file import *


class GenerateTest:
    def __init__(self):
        self.sr = SimpleRead()

    def generate_from_ls(self, test_filename) -> list:
        return self.sr.read_each_line(test_filename)

    @staticmethod
    def generate_from_composition(prefix: str, infix: list, suffix=""):
        return [prefix + inf + suffix for inf in infix]
