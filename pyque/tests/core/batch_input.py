#!/usr/bin/env python3
# created on Jan 8, 2018 at 19:20 by Qi Zhang

from os.path import join, pardir
from unittest import TestCase

from pyque.core.batch_input import BatchInput
from pyque.lexer.batch import BatchTemplateLexer


class BatchInputTester(TestCase):
    def setUp(self):
        self.infile = join(pardir, 'data/sample_batch_template')
        self.template_dict = BatchTemplateLexer(infile=self.infile).to_dict()
        self.bi = BatchInput('slurm')

    def test_shebang(self):
        self.bi.shebang = self.template_dict['shebang']
        self.assertEqual(self.bi.shebang, '#!/bin/sh')

    def test_modules(self):
        for module in self.template_dict['modules']:
            self.bi.add_modules(module)
        print(self.bi.modules)

    def test_build_scheduler(self):
        self.bi.scheduler.account = self.template_dict['account']
        print(self.bi.scheduler.account)
