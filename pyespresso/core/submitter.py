#!/usr/bin/env python3
"""
:mod:`pw` -- title
========================================

.. module dicts
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import shlex
import subprocess


class Submitter:
    def __init__(self, batch_file: str):
        self.batch_file = batch_file
        self.process = None
        self.args = None

    def submit(self):
        if not self.args:
            args = ['sbatch', self.batch_file]
        else:
            args = ['sbatch'] + shlex.split(self.args) + [self.batch_file]
        with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
            self.process = proc
            return self.get_job_id()

    def get_job_id(self) -> str:
        return str(self.process.stdout.read().split()[-1])
