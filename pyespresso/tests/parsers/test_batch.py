#!/usr/bin/env python3

from os.path import join, pardir
from unittest import TestCase

from pyespresso.lexer.batch import BatchTemplateLexer


class BatchTemplateParserTester(TestCase):
    def setUp(self):
        self.infile = join(pardir, 'data/sample_batch_template')
        self.btp = BatchTemplateLexer(infile=self.infile)

    def test_to_dict(self):
        print(self.btp.to_dict())
