#!/usr/bin/env python3
# created at Dec 22, 2017 10:48 PM by Qi Zhang

from readers.phonon import PHononInputParser
from basics.phonon import PHononStandaradInput


class PHononInputBuilder:
    def __init__(self, in_file):
        self.parser = PHononInputParser(in_file)
        self.input_obj = PHononStandaradInput()

    def build_INPUTPH(self):
        self.input_obj.INPUTPH_namelist = self.parser.parse_INPUTPH_namelist()

    def build_title(self):
        self.input_obj.__title__ = self.parser.parse_title()

    def build_single_q_point(self):
        self.input_obj.single_q_point = self.parser.parse_single_q_point()

    def build_q_points(self):
        self.input_obj.q_points = self.parser.parse_q_points()

    def build_all(self):
        self.build_INPUTPH()
        # self.build_title()
        self.build_single_q_point()
        self.build_q_points()
