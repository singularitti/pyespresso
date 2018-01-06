#!/usr/bin/env python3
# created on Jan 5, 2018 at 23:29 by Qi Zhang

from unittest import TestCase
from meta.text import TextStream


class StreamTester(TestCase):
    def setUp(self):
        self.infile = 'cell_test.py'
        with open(self.infile, 'r') as f:
            self.instream = repr(f.read())
        self.stream = TextStream(instr=self.instream)

    def test_print(self):
        for line in self.stream.stream_generator():
            print(line)

    def test_to_string_io(self):
        print(next(self.stream.to_string_io))
