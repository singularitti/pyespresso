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

import smtplib
import subprocess
import time
from typing import *

from lazy_property import LazyWritableProperty

from pyespresso.core.email_sender import EmailSender


class SlurmSystemMonitor:
    def __init__(self, job_id: Union[str, int]):
        self.job_id = str(job_id)

    def cancel(self):
        subprocess.call(['scancel', self.job_id])

    @staticmethod
    def share():
        subprocess.call(['sshare'])

    def control(self):
        subprocess.call(['scontrol', 'show', 'job', self.job_id])

    @LazyWritableProperty
    def user_name(self) -> str:
        pass

    def queue_jobs_by_user_name(self):
        if not self.user_name:
            raise AttributeError("'user_name' attribute is not set!")
        subprocess.call(['squeue', '-u', self.user_name])

    @LazyWritableProperty
    def account(self):
        pass

    def queue_jobs_by_account(self):
        if not self.account:
            raise AttributeError("'account' attribute is not set!")
        subprocess.call(['squeue', '-A', self.account])

    def queue_job_by_id(self):
        subprocess.call(['squeue', '--job', self.job_id])

    @LazyWritableProperty
    def process(self):
        pass

    def is_process_alive(self):
        job_status = subprocess.Popen(['squeue', '-j', self.job_id], stdout=subprocess.PIPE)
        queue = job_status.stdout.read().split()
        if self.job_id in queue:
            print("Congratulations, the variable cell relaxation is done!")
            return True
        return False

    def queue_every_n_seconds(self, n: float):
        if self.is_process_alive():
            print(self.queue_job_by_id())
        else:
            print('Process {0} is end!'.format(self.job_id))
        time.sleep(n)

    @LazyWritableProperty
    def user_email(self):
        pass

    @LazyWritableProperty
    def user_email_password(self):
        pass

    @LazyWritableProperty
    def job_name(self):
        pass

    def email_me_after_done(self, n: float):
        if not self.is_process_alive():
            email_sender = EmailSender(self.user_email, self.user_email, 'Job {0} is done!'.format(self.job_id))
            email_sender.server = smtplib.SMTP('smtp.gmail.com', 587)
            email_sender.message = "Your job with name {0} and job ID {1} is done!".format(self.job_name, self.job_id)
            email_sender.from_address_password = self.user_email_password
            email_sender.send()
        else:
            time.sleep(n)
