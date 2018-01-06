#!/usr/bin/env python3
# created on Jan 6, 2018 at 02:18 by Qi Zhang

from submitters.scf_generator import PWscfStandardInputCreator
from unittest import TestCase
from beeprint import pp

instream = \
    """
    pseudopotential directory: /rigel/home/qz2280/pseudo
    prefix: hcpFe
    fictitious cell mass: 0.00010
    number of atoms: 2
    number of atomic types: 1
    energy cutoff: 90
    smearing: fd
    occupations: smearing
    degauss: 0.0018874190668942044
    scratch folder: /rigel/mphys/users/qz2280/vcrelax/test31
    """


class PWscfStandardInputCreatorTester(TestCase):
    def setUp(self):
        self.creator = PWscfStandardInputCreator(instream=instream)

    def test_read_fixed_input(self):
        d = self.creator.build_parameters()
        pp(d)
