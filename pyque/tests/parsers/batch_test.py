#!/usr/bin/env python3
# created on Jan 8, 2018 at 18:56 by Qi Zhang

from os.path import join, pardir
from unittest import TestCase

from pyque.lexer.batch import BatchTemplateParser


class BatchTemplateParserTester(TestCase):
    def setUp(self):
        self.infile = join(pardir, 'data/sample_batch_template')
        self.btp = BatchTemplateParser(infile=self.infile)

    def test_to_dict(self):
        print(self.btp.to_dict())
