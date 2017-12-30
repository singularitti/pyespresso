#!/usr/bin/env python3
# created at Dec 30, 2017 1:14 PM by Qi Zhang


# Referenced from
# [here](http://nbviewer.jupyter.org/urls/gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb).
class Descriptor(object):
    def __init__(self):
        # notice we aren't setting the label here
        self.label = None

    def __get__(self, instance, owner):
        return instance.__dict__.get(self.label, None)

    def __set__(self, instance, value):
        instance.__dict__[self.label] = value


class DescriptorOwner(type):
    def __new__(cls, name, bases, attrs):
        # find all descriptors, auto-set their labels
        for n, v in attrs.items():
            if isinstance(v, Descriptor):
                v.label = n
        return super(DescriptorOwner, cls).__new__(cls, name, bases, attrs)
