#!/usr/bin/env python3
# created at Dec 24, 2017 7:34 PM by Qi Zhang

import os
from typing import Iterable, Union, Set, List

from lazy_property import *

from miscellaneous.sets import add_elements_to_set, remove_elements_from_set
from submitters.scheduler import available_schedulers, SchedulerSystem


def _is_any_not_string(iterable: Iterable) -> bool:
    """
    If any element of an iterable is not a string, return `True`.

    :param iterable: Can be a set, a tuple, a list, etc.
    :return: Whether any element of an iterable is not a string.
    """
    return any(type(_) != str for _ in iterable)


# ========================================= These are some useful classes. =========================================
class Comment:
    """
    Insert comment between lines.
    """

    def __init__(self, raw_comment: str, line_length: int = 118, comment_character: str = '#'):
        """

        :param raw_comment: Comment string given by user.
        :param line_length: Each line of the comment has a fixed length of characters.
        :param comment_character:
        """
        self.raw_comment = raw_comment
        self.fixed_length = line_length
        self.comment_character = comment_character

    @staticmethod
    def chunk_string(string, length):
        return (string[0 + i:length + i] for i in range(0, len(string), length))

    def __str__(self):
        lines = (_.strip() for _ in self.raw_comment.splitlines())
        comment = []
        for line in lines:
            for chunk in self.chunk_string(line, self.fixed_length):
                comment.append(self.comment_character + ' ' + chunk)
        return "\n".join(comment)

    __repr__ = __str__


_comment = Comment  # Alias


# ========================================= The following are core classes. =========================================
class JobInput:
    def __init__(self, scheduler_name: str, directive_style: str = 'short'):
        """

        :param scheduler_name:
        """
        self.__name__ = 'JobHeader'
        self.__directive_style: str = directive_style
        self.scheduler: SchedulerSystem = available_schedulers[scheduler_name]()
        self.__modules: Set = set()
        self.__commands: List[str] = []

    @LazyProperty
    def scheduler_name(self) -> str:
        """

        :return: 'Slurm system' or 'Cray system'
        """
        return self.scheduler.__name__

    @property
    def directive_style(self):
        return self.__directive_style

    @directive_style.setter
    def directive_style(self, new_style: str):
        if new_style.lower() not in {'short', 'long'}:
            raise ValueError()
        else:
            self.__directive_style = new_style

    @LazyWritableProperty
    def shebang(self) -> str:
        pass

    @property
    def modules(self) -> Set[str]:
        return self.__modules

    def add_modules(self, *args: Union[str, Iterable[str]]) -> None:
        if _is_any_not_string(args):
            raise TypeError('Modules added should all be strings! Check your type!')
        if not args:  # If nothing is provided
            raise ValueError('No modules are added!')
        self.__modules = add_elements_to_set(self.__modules, args)
        if len(args) == 1:
            print('Module {0} is added!'.format(args[0]))
        else:  # len(args) > 1
            print("Modules {0} are added!".format(', '.join(args)))

    def remove_modules(self, *args) -> None:
        if _is_any_not_string(args):
            raise TypeError('Modules removed should all be strings! Check your type!')
        self.__modules = remove_elements_from_set(self.__modules, args)

    def collect_short_directives(self):
        directives = dict()
        for obj in type(self.scheduler).__dict__.values():
            if hasattr(obj, 'short_directive'):
                print((obj))
                directives.update({type(obj).__name__: (obj.short_directive, obj.__get__(self.scheduler, None))})
        return directives

    def collect_long_directives(self):
        directives = dict()
        for obj in dir(type(self.scheduler)):
            if hasattr(obj, 'long_directive'):
                directives.update({type(obj).__name__: (obj.long_directive, self.scheduler.__dict__[type(obj).__name__])})
        return directives

    @staticmethod
    def add_comment(raw_comment: str) -> str:
        return str(_comment(raw_comment))

    @property
    def commands(self):
        return self.__commands

    def add_command(self, new_command):
        self.__commands.append(new_command)

    def remove_command(self, command):
        self.__commands.remove(command)

    def write_to_file(self, output_name: str, output_path: str = '') -> None:
        """

        :param output_name: A path redirects to the output file you want.
        :param output_path:
        :return:
        """
        if not output_path:
            file_path = os.path.join(os.getcwd(), output_name)
        else:
            file_path = os.path.join(output_path, output_name)

        scheduler = self.scheduler
        directive_prefix = scheduler.directive_prefix
        with open(file_path, 'w') as f:
            f.write("{0}\n".format(self.shebang))
            if self.directive_style == 'short':
                print(self.collect_short_directives().items())
                for directive_name, directive in self.collect_short_directives().items():
                    f.write("{0} {1}={2}\n".format(directive_prefix, directive[0], directive[1]))
            else:
                for directive_name, directive in self.collect_long_directives().items():
                    f.write("{0} {1}={2}".format(directive_prefix, directive[0], directive[1]))
            f.write("\n")
            for module in self.modules:
                f.write("module load {0}\n".format(module))
            f.write("\n")
            f.write("\n".join(self.commands))

        print("A job input has been written to file {0}!".format(file_path))
