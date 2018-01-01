#!/usr/bin/env python3
# created at Dec 30, 2017 1:14 PM by Qi Zhang

from weakref import WeakKeyDictionary


# Referenced from
# [here](http://nbviewer.jupyter.org/urls/gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb).
class Descriptor(object):
    """
    A "strong" descriptor.
    """

    def __init__(self, default=None):
        # notice we aren't setting the label here
        self.label = None
        self.default = default

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return instance.__dict__.get(self.label, self.default)

    def __set__(self, instance, value):
        instance.__dict__[self.label] = value


class DescriptorOwner(type):
    def __new__(cls, name, bases, attrs):
        # find all descriptors, auto-set their labels
        for n, v in attrs.items():
            if isinstance(v, Descriptor):
                v.label = n
        return super(DescriptorOwner, cls).__new__(cls, name, bases, attrs)


class WeakDescriptor(object):
    def __init__(self, default=None):
        # notice we aren't setting the label here
        self.data = WeakKeyDictionary()
        self.default = default

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return self.data.get(instance, self.default)

    def __set__(self, instance, value):
        self.data[instance] = value
