#!/usr/bin/env python3
"""
:mod:`pw` -- title
========================================

.. module dicts
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

import subprocess
from lazy_property import LazyWritableProperty
import shlex


class Submitter:
    def __init__(self, batch_file: str):
        self.batch_file = batch_file
        self.process = None

    @LazyWritableProperty
    def args(self):
        pass

    def submit(self):
        if not self.args:
            args = ['sbatch', self.batch_file]
        else:
            args = ['sbatch'] + shlex.split(self.args) + [self.batch_file]
        with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
            self.process = proc
            return self.get_job_id()

    def get_job_id(self):
        return int(self.process.stdout.read().split()[-1])

    # def _distribute_task(self) -> int:
    #     """
    #     Distribute calculations on different pressures evenly to all processors you have.
    #
    #     :return: number of cores per pressure
    #     """
    #     print('You have {0} pressures and {1} cores.'.format(self.pressures_num, self.cores_num * self.nodes_num))
    #     div = int(self.nodes_num * self.cores_num / self.pressures_num)
    #
    #     if (self.nodes_num * self.cores_num) % self.pressures_num == 0:
    #         print('Therefore, I will consider {0} cores per pressure.'.format(div))
    #     else:
    #         raise TypeError('Number of cores is not divided by number of pressures!')
    #
    #     return div
