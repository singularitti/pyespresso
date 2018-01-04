#!/usr/bin/env python3
# created at Nov 25, 2017 8:05 PM by Qi Zhang
from unittest import TestCase

from beeprint import pp

from data_models.phonon_builder import PHononInputBuilder
from readers.phonon import *


class TestPHononInputReader(TestCase):
    def setUp(self):
        self.phip = PHononInputParser('ph.in')

    def test_parse_INPUTPH_namelist(self):
        pp(self.phip.parse_INPUTPH_namelist())


class TestPHononInputBuilder(TestCase):
    def setUp(self):
        self.phib = PHononInputBuilder('ph.in')

    def test_build_all(self):
        pp(self.phib.input_obj.single_q_point)
        self.phib.build_all()
