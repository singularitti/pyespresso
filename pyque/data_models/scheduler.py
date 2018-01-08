#!/usr/bin/env python3
# created at Dec 29, 2017 9:17 PM by Qi Zhang

import warnings
from time import strptime
from typing import Union

from lazy_property import LazyWritableProperty, LazyProperty

from default_configuerations.scheduler import DEFAULT_SLURM_CONFIG
from meta.descriptors import LabeledDescriptor, MetaDescriptorOwner


class _NodesNumber(LabeledDescriptor):
    def __set__(self, instance, new_nodes_number: int):
        if type(new_nodes_number) is not int:
            raise TypeError('The number of nodes is not an integer!')
        if new_nodes_number < 0:
            raise ValueError('Wrong number of nodes is given!')
        # You need not to put it in an `else` clause because the program will end if it meets an error.
        instance.__dict__[self.label] = new_nodes_number


class _Time(LabeledDescriptor):
    def __set__(self, instance, new_time: str):
        if instance is None:
            raise AttributeError
        if type(new_time) is not str:
            raise TypeError('The time you set is not a string!')
        hour: int = strptime(new_time, '%H:%M:%S').tm_hour
        if hour > 120:
            warnings.warn('The maximum time allowed is five days!')
        instance.__dict__[self.label] = new_time


def _slurm_system_config(item: str):
    if True:
        return DEFAULT_SLURM_CONFIG[item]


class SlurmSystem(metaclass=MetaDescriptorOwner):
    def __init__(self):
        self.__name__ = 'Slurm system'

    # Class level attribute
    nodes_number = _NodesNumber(_slurm_system_config('nodes_number'))
    nodes_number.long_directive = '--nodes'
    nodes_number.short_directive = '-N'
    N = nodes_number  # Attribute alias

    time = _Time(_slurm_system_config('time'))
    time.long_directive = '--time'
    time.short_directive = '-t'
    t = time  # Attribute alias

    @LazyProperty
    def directive_prefix(self):
        """
        A read-only attribute.

        :return:
        """
        return '#SBATCH'

    @LazyWritableProperty
    def account(self) -> str:
        pass

    account.long_directive = '--account'
    account.short_directive = '-A'

    A = account  # Attribute alias

    @LazyWritableProperty
    def job_name(self) -> str:
        pass

    job_name.long_directive = '--job-name'
    job_name.short_directive = '-J'

    J = job_name  # Attribute alias

    @LazyWritableProperty
    def mail_type(self) -> str:
        pass

    # mail_type.long_directive = '--mail-type'

    @LazyWritableProperty
    def mail_user(self) -> str:
        pass

    # mail_user.long_directive = '--mail-user'

    @LazyWritableProperty
    def mem(self) -> int:
        pass

    # mem.long_directive = '--mem'

    @LazyWritableProperty
    def tasks_per_node(self) -> int:
        pass

    @LazyWritableProperty
    def cpus_per_task(self) -> int:
        pass

    # cpus_per_task.long_directive = '--cpus-per-task'

    @LazyProperty
    def exclusive(self):
        pass


available_schedulers = {
    'slurm': SlurmSystem
}

SchedulerSystem = Union[SlurmSystem]
