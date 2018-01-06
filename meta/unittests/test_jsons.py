#!/usr/bin/env python3
# created on Jan 6, 2018 at 08:48 by Qi Zhang

from unittest import TestCase

from meta.jsons import JsonStr

instr = \
    """
    {
      "&control": {
        "calculation": "scf",
        "prefix": "silicon",
        "restart_mode": "from_scratch",
        "tstress": true,
        "tprnfor": true,
        "pseudo_dir": "./",
        "outdir": "./"
      },
      "&system": {
        "ibrav": 2,
        "celldm(1)": 10.2,
        "nat": 2,
        "ntyp": 1,
        "ecutwfc": 18.0
      },
      "&electrons": {},
      "ATOMIC_SPECIES": {
        "Si": [28.086, "Si.vbc.UPF"]
      },
      "ATOMIC_POSITIONS": {
        "Si": [
          [0.00, 0.00, 0.00],
          [0.25, 0.25, 0.25]
        ]
      },
      "K_POINTS": {
        "nk": [4, 4, 4],
        "sk": [1, 1, 1]
      }
    }
    """


class JsonStrTester(TestCase):
    def setUp(self):
        self.jsonstr = JsonStr(instr)
        self.json2 = JsonStr(infile='si.scf.json')

    def test_content(self):
        print(self.jsonstr.content())
        print(self.json2.content())
