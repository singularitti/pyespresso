#!/usr/bin/env python3
# created at Dec 29, 2017 9:17 PM by Qi Zhang

from time import strptime
from typing import Union

from lazy_property import LazyWritableProperty

from configue.scheduler import default_Slurm_config
from miscellaneous.descriptors import Descriptor, DescriptorOwner


class _NodesNumber(Descriptor):
    def __init__(self, default_nodes_number: int):
        super().__init__(default_nodes_number)
        self.__name__ = 'nodes'
        self.has_short_directive = True
        self.long_directive_suffix = '--nodes'
        self.short_directive_suffix = '-N'

    def __set__(self, instance, new_nodes_number: int):
        if type(new_nodes_number) is not int:
            raise TypeError('The number of nodes is not an integer!')
        if new_nodes_number < 0:
            raise ValueError('Wrong number of nodes is given!')
        # You need not to put it in an `else` clause because the program will end if it meets an error.
        instance.__dict__[self.label] = new_nodes_number


class _Time(Descriptor):
    def __init__(self, default_time: str):
        super().__init__(default_time)
        self.__name__ = 'time'
        self.has_short_directive = True
        self.long_directive_suffix = '--time'
        self.short_directive_suffix = '-t'

    def __set__(self, instance, new_time: str):
        if type(new_time) is not str:
            raise TypeError('The time you set is not a string!')
        hour: int = strptime(new_time, '%H:%M:%S').tm_hour
        if hour > 120:
            print('The maximum time allowed is five days!')
        instance.__dict__[self.label] = new_time


def _slurm_system_config(item: str):
    if True:
        return default_Slurm_config[item]


class SlurmSystem(metaclass=DescriptorOwner):
    nodes_number = _NodesNumber(_slurm_system_config('nodes_number'))
    N = nodes_number  # Attribute alias

    time = _Time(_slurm_system_config('time'))
    t = time  # Attribute alias

    def __init__(self):
        self.__name__ = 'Slurm system'

    @property
    def directive_prefix(self):
        """
        A read-only attribute.

        :return:
        """
        return '#SBATCH'

    @LazyWritableProperty
    def account(self) -> str:
        pass

    A = account  # Attribute alias

    @LazyWritableProperty
    def job_name(self) -> str:
        pass

    J = job_name  # Attribute alias

    @LazyWritableProperty
    def mail_type(self) -> str:
        pass

    @LazyWritableProperty
    def mail_user(self) -> str:
        pass

    @LazyWritableProperty
    def mem(self) -> int:
        pass

    @LazyWritableProperty
    def tasks_per_node(self) -> int:
        pass


SlurmSystem.account.__name__ = 'account'
SlurmSystem.account.long_directive_suffix = '--account'
SlurmSystem.account.short_directive_suffix = '-A'

SlurmSystem.job_name.__name__ = 'job-name'
SlurmSystem.job_name.long_directive_suffix = '--job-name'
SlurmSystem.job_name.short_directive_suffix = '-J'

available_schedulers = {
    'slurm': SlurmSystem
}

SchedulerSystem = Union[SlurmSystem]
