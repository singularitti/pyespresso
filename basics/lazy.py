#!/usr/bin/env python3
# created at Dec 2, 2017 5:25 PM by Qi Zhang

from functools import update_wrapper


class _Missing(object):

    def __repr__(self):
        return 'no value'

    def __reduce__(self):
        return '_missing'


_missing = _Missing()


class CachedProperty(property):
    """
    A decorator that converts a function into a lazy property.  The
    function wrapped is called the first time to retrieve the result
    and then that calculated result is used the next time you access
    the value::

        class Foo(object):
            @CachedProperty
            def foo(self):
                # calculate something important here
                return 42

    The class has to have a `__dict__` in order for this property to
    work.

    Referenced from [here](https://github.com/pallets/werkzeug/blob/master/werkzeug/utils.py).
    """

    # implementation detail: this property is implemented as non-data
    # descriptor.  non-data descriptors are only invoked if there is
    # no entry with the same name in the instance's __dict__.
    # this allows us to completely get rid of the access function call
    # overhead.  If one chooses to invoke __get__ by hand the property
    # will still work as expected because the lookup logic is replicated
    # in __get__ for manual invocation.

    def __init__(self, func, name=None, doc=None):
        self.__name__ = name or func.__name__
        self.__module__ = func.__module__
        self.__doc__ = doc or func.__doc__
        self.func = func

    def __set__(self, obj, value):
        obj.__dict__[self.__name__] = value

    def __get__(self, obj, type=None):
        if obj is None:
            return self
        value = obj.__dict__.get(self.__name__, _missing)
        if value is _missing:
            value = self.func(obj)
            obj.__dict__[self.__name__] = value
        return value


class LazyProperty(property):
    def __init__(self, method, fget=None, fset=None, fdel=None, doc=None):

        self.method = method
        self.cache_name = "_{}".format(self.method.__name__)

        doc = doc or method.__doc__
        super(LazyProperty, self).__init__(fget=fget, fset=fset, fdel=fdel, doc=doc)

        update_wrapper(self, method)

    def __get__(self, instance, owner):
        """

        :param instance: Either an instance `b` of a class (for `b.x`) or `None` (for `B.x`).
        :param owner: A class
        :return:
        """

        if instance is None:
            return self  # `self` is `B.x`, `B` is the class, `x` is the descriptor (the `LazyProperty`)

        if hasattr(instance, self.cache_name):  # instance is `b`, does it have like `_account` attribute
            result = getattr(instance, self.cache_name)
        else:
            if self.fget is not None:
                result = self.fget(instance)
            else:
                result = self.method(instance)

            setattr(instance, self.cache_name, result)

        return result


class LazyWritableProperty(LazyProperty):
    def __set__(self, instance, value):

        if instance is None:  # D.m.__set__(None, 1) will raise Error
            raise AttributeError

        if self.fset is None:
            setattr(instance, self.cache_name, value)
        else:
            self.fset(instance, value)
