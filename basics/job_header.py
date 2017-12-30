#!/usr/bin/env python3
# created at Dec 24, 2017 7:34 PM by Qi Zhang

from typing import Iterable, Union, Set

from miscellaneous.sets import add_elements_to_set, remove_elements_from_set


def _is_any_not_string(iterable: Iterable) -> bool:
    """
    If any element of an iterable is not a string, return `True`.

    :param iterable: Can be a set, a tuple, a list, etc.
    :return: Whether any element of an iterable is not a string.
    """
    return any(type(_) != str for _ in iterable)


class JobHeader:
    def __init__(self):
        self.__name__ = 'JobHeader'
        self._account: str = ''
        self._job_name: str = ''
        self._scheduler: str = ''
        self._modules: Set = set()
        self._time: str = '00:00:00'
        self._nodes_num = 0
        self._processors_num: int = 0
        self._tasks_per_node = 0
        self._shebang = ''

    @property
    def account(self) -> str:
        if self._account == '':
            print("No account name is given! Use default value ''!")
        return self._account

    @account.setter
    def account(self, new_account: str) -> None:
        self._account = new_account

    @property
    def job_name(self):
        if self._job_name == '':
            print("No job name is given! Use default value ''!")
        return self._job_name

    @job_name.setter
    def job_name(self, new_job_name: str):
        self._job_name = new_job_name

    @property
    def scheduler(self) -> str:
        return self._scheduler

    @scheduler.setter
    def scheduler(self, new_scheduler: str):
        if new_scheduler.upper() not in {'SLURM', 'CRAY'}:
            raise ValueError('Unknown scheduler {0} is provided!'.format(new_scheduler))
        self._scheduler = new_scheduler

    @property
    def modules(self) -> Set[str]:
        return self._modules

    @modules.setter
    def modules(self, new_modules: Set[str]) -> None:
        if _is_any_not_string(new_modules):
            raise TypeError('Modules should all be strings! Check your type!')
        self._modules = new_modules

    def add_modules(self, *args: Union[str, Iterable[str]]) -> None:
        if _is_any_not_string(args):
            raise TypeError('Modules added should all be strings! Check your type!')
        if not args:  # If nothing is provided
            raise ValueError('No modules are added!')
        self._modules = add_elements_to_set(self._modules, args)
        if len(args) == 1:
            print('Module {0} is added!'.format(args[0]))
        else:  # len(args) > 1
            print("Modules {0} are added!".format(', '.join(args)))

    def remove_modules(self, *args) -> None:
        if _is_any_not_string(args):
            raise TypeError('Modules removed should all be strings! Check your type!')
        self._modules = remove_elements_from_set(self._modules, args)

    @property
    def time(self) -> str:
        return self._time

    @time.setter
    def time(self, new_time: str) -> None:
        self._time = new_time

    @property
    def processors_num(self) -> int:
        return self._processors_num

    @processors_num.setter
    def processors_num(self, new_processors_num: int) -> None:
        if type(new_processors_num) == int:
            raise TypeError('The number of processors should be an integer!')
        self._processors_num = new_processors_num

    @property
    def nodes_num(self) -> int:
        return self._nodes_num

    @nodes_num.setter
    def nodes_num(self, new_nodes_num: int) -> None:
        if type(new_nodes_num) == int:
            raise TypeError('The number of nodes should be an integer!')
        self._nodes_num = new_nodes_num

    @property
    def tasks_per_node(self) -> int:
        return self._tasks_per_node

    @tasks_per_node.setter
    def tasks_per_node(self, new_tasks_per_node):
        self._tasks_per_node = new_tasks_per_node

    @property
    def shebang(self) -> str:
        return self._shebang

    @shebang.setter
    def shebang(self, new_shebang: str) -> None:
        self._shebang = new_shebang

    def write_to_file(self, output_file: str) -> None:
        """


        :param output_file: A path redirects to the output file you want.
        :return:
        """
        with open(output_file, 'w') as f:
            f.write("{0}\n".format(self.shebang))
