#!/usr/bin/env python3
# created at Dec 24, 2017 7:34 PM by Qi Zhang

from typing import Iterable, Union, Set

from lazy_property import LazyWritableProperty

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
class JobHeader:
    def __init__(self, scheduler_name: str, directive_style: str = 'short'):
        """

        :param scheduler_name:
        """
        self.__name__ = 'JobHeader'
        self.scheduler_name: str = scheduler_name
        self.directive_style: str = directive_style
        self.__scheduler: SchedulerSystem = available_schedulers[scheduler_name]
        self.__modules: Set = set()

    @LazyWritableProperty
    def shebang(self) -> str:
        pass

    @property
    def modules(self) -> Set[str]:
        return self.__modules

    @modules.setter
    def modules(self, new_modules: Set[str]) -> None:
        if _is_any_not_string(new_modules):
            raise TypeError('Modules should all be strings! Check your type!')
        self.__modules = new_modules

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

    def collect_directives(self):
        directives = dict()
        for obj in vars(self.__scheduler):
            if obj.has_short_directive:
                if self.directive_style == 'short':
                    directives.update({obj: obj.short_directive})
                else:  # self.directive_style == 'long'
                    directives.update({obj: obj.long_directive})

        return directives

    @staticmethod
    def add_comment(raw_comment: str) -> str:
        return str(_comment(raw_comment))

    @LazyWritableProperty
    def commands(self):
        pass

    def write_to_file(self, output_file: str) -> None:
        """


        :param output_file: A path redirects to the output file you want.
        :return:
        """
        directive_prefix = self.__scheduler.directive_prefix
        with open(output_file, 'w') as f:
            f.write("{0}\n".format(self.shebang))
            for directive in self.collect_directives():
                f.write(
                    "{0} {1}={2}".format(directive_prefix, directive, self.__scheduler.__dict__[directive.__name__]))
            for module in self.modules:
                f.write("module load {0}\n".format(module))
