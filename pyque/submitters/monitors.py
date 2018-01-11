#!/usr/bin/env python3
"""
:mod:`monitors` -- 
========================================

.. module dicts
   :platform: Unix, Windows, Mac, Linux
   :synopsis: We use the ``subprocess`` module instead of ``os.system`` because the former one does shell escaping
   and is therefore much safer.
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

from subprocess import call
from typing import *

from lazy_property import LazyWritableProperty


class SlurmSystemMonitor:
    def __init__(self, job_id: Union[str, int]):
        self.job_id = str(job_id)

    def cancel(self):
        call(['scancel', self.job_id])

    @staticmethod
    def share():
        call(['sshare'])

    def control(self):
        call(['scontrol', 'show', 'job', self.job_id])

    @LazyWritableProperty
    def user_name(self) -> str:
        pass

    def queue_jobs_by_user_name(self):
        if not self.user_name:
            raise AttributeError("'user_name' attribute is not set!")
        call(['squeue', '-u', self.user_name])

    @LazyWritableProperty
    def account(self):
        pass

    def queue_jobs_by_account(self):
        if not self.account:
            raise AttributeError("'account' attribute is not set!")
        call(['squeue', '-A', self.account])

    def queue_job_by_id(self):
        call(['squeue', '--job', self.job_id])
