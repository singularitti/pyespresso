#!/usr/bin/env python3
# created on Jan 6, 2018 at 02:18 by Qi Zhang

from unittest import TestCase

from beeprint import pp

from submitters.scf_generator import PWscfFixedFormParser

instream = \
    """
    pseudopotential directory: ./
    prefix: silicon
    fictitious cell mass: 0.00010
    number of atoms: 2
    number of atomic types: 1
    energy cutoff: 18.0
    smearing: fd
    occupations: smearing
    degauss: 0.0018874190668942044
    scratch folder: ./
    """


class PWscfStandardInputCreatorTester(TestCase):
    def setUp(self):
        self.creator = PWscfFixedFormParser(instream=instream)

    def test_raw_parameters(self):
        pp(self.creator.raw_parameters)

    def test_build_parameters(self):
        pp(self.creator._build_parameters())
