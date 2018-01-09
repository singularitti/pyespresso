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
import shlex
from lazy_property import LazyWritableProperty


class Submitter:
    # def submit(self):
    #     with subprocess.Popen(self.args, stdout=subprocess.PIPE) as proc:
    #         job_id = int(proc.stdout.read().split()[-1])
    #         print('Job submitted. The job ID is {0}.'.format(job_id))
    #         print('Waiting for calculation. This can take a while, you can grab a coffee!')
    #         self.job_id = job_id
    #         return job_id

    def _distribute_task(self) -> int:
        """
        Distribute calculations on different pressures evenly to all processors you have.

        :return: number of cores per pressure
        """
        print('You have {0} pressures and {1} cores.'.format(self.pressures_num, self.cores_num * self.nodes_num))
        div = int(self.nodes_num * self.cores_num / self.pressures_num)

        if (self.nodes_num * self.cores_num) % self.pressures_num == 0:
            print('Therefore, I will consider {0} cores per pressure.'.format(div))
        else:
            raise TypeError('Number of cores is not divided by number of pressures!')

        return div
