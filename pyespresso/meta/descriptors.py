#!/usr/bin/env python3
# created at Dec 30, 2017 1:14 PM by Qi Zhang
"""
Referenced from
`here <http://nbviewer.jupyter.org/urls/gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb/>`_.
"""

from weakref import WeakKeyDictionary


class LabeledDescriptor:
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


class DescriptorOwnerMeta(type):
    def __new__(cls, name, bases, attrs):
        # find all descriptors, auto-set their labels
        for n, v in attrs.items():
            if isinstance(v, LabeledDescriptor):
                v.label = n
        return super(DescriptorOwnerMeta, cls).__new__(cls, name, bases, attrs)


class WeakDescriptor:
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
