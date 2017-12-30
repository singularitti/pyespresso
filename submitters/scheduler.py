#!/usr/bin/env python3
# created at Dec 29, 2017 9:17 PM by Qi Zhang

from time import strptime

from lazy_property import LazyProperty, LazyWritableProperty

from configue.scheduler import default_Slurm_config
from miscellaneous.descriptors import Descriptor, DescriptorOwner


class _NodesNumber(Descriptor):
    def __init__(self, default_nodes_number: int):
        super().__init__()
        self.default_nodes_number = default_nodes_number

    def __get__(self, instance, owner):
        return instance.__dict__.get(self.label, self.default_nodes_number)

    def __set__(self, instance, new_nodes_number: int):
        if type(new_nodes_number) is not int:
            raise TypeError('The number of nodes is not an integer!')
        if new_nodes_number < 0:
            raise ValueError('Wrong number of nodes is given!')
        else:
            instance.__dict__[self.label] = new_nodes_number


class _Time(Descriptor):
    def __init__(self, default_time: str):
        super().__init__()
        self.default_time = default_time

    def __get__(self, instance, owner):
        return instance.__dict__.get(self.label, self.default_time)

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
    time = _Time(_slurm_system_config('time'))

    def __init__(self):
        self.__name__ = 'Slurm system'

    @LazyProperty
    def directive_prefix(self):
        return '#SBATCH'

    @LazyWritableProperty
    def account(self) -> str:
        pass

    @LazyWritableProperty
    def job_name(self) -> str:
        pass

    @LazyWritableProperty
    def mail_type(self) -> str:
        pass

    @LazyWritableProperty
    def mail_user(self) -> str:
        pass

    @LazyWritableProperty
    def mem(self) -> int:
        pass
