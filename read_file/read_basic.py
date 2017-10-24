#!/usr/bin/env python3
# created at Oct 20, 2017 6:13 PM by Qi Zhang

import re
from typing import *


class SimpleReader:
    def __init__(self, in_file: str):
        """
        For each file, you need to generate a `SimpleRead` class to read them.

        :param in_file:
        """
        self.in_file = in_file

    def _match_only_once(self, pattern: str, wrapper: Callable):
        """
        This method matches a pattern which exists only once in a file.

        :param pattern: a regular expression that you want to match
        :param wrapper: a wrapper function which determines the returned type of value
        :return: Determined by the `wrapper`, the value you want to grep out from the file.
        """
        with open(self.in_file, 'r') as f:
            match = re.findall(pattern, f.read())
        return wrapper(match[0])

    def read_each_line(self) -> List[str]:
        """
        This method reads each line simply from a file, and discard the '\n' character.

        :return: A list contains all the lines of the file.
        """
        line_list = []
        with open(self.in_file, 'r') as f:
            for line in f:
                line_list.append(re.split('\n', line)[0])
        return line_list

    def read_two_columns(self) -> Tuple[List[str], List[str]]:
        """
        This method reads 2 columns from a file.

        :return: two columns of the file
        """
        col1_list = []
        col2_list = []
        with open(self.in_file, 'r') as f:
            for line in f:
                if not line.split():
                    f.readline()
                else:
                    sp = line.split()
                    col1_list.append(sp[0])
                    col2_list.append(sp[1])
        return col1_list, col2_list

    def read_one_column_as_keys(self, col_index: int, wrapper: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param col_index: the index of the column that you want to make it as keys.
        :param wrapper: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        key_list = []
        value_list = []
        # Add utf-8 support because we may use special characters.
        with open(self.in_file, 'r', encoding='utf-8') as f:
            for line in f:
                sp = line.split()
                key_list.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                value_list.append(wrapper(sp))
        return dict(zip(key_list, value_list))

    def _read_reciprocal_points(self) -> Dict[str, List[float]]:
        """
        Suppose you have a file like this:
            A	0.0000000000	0.0000000000	0.5000000000
            Γ   0.0000000000	0.0000000000	0.0000000000
            H	0.3333333333	0.3333333333	0.5000000000
            H2	0.3333333333	0.3333333333   -0.5000000000
            K	0.3333333333	0.3333333333	0.0000000000
            L	0.5000000000	0.0000000000	0.5000000000
            M	0.5000000000	0.0000000000	0.0000000000
        These are the k-points you want to track through.
        This method reads through those names and numbers, and set each name as a key, each 3 k-coordinates as
        its value, forms a dictionary.

        :return: a dictionary
        """
        return self.read_one_column_as_keys(0, lambda x: list(map(float, x)))


def _str_list_to(inp: List[str], to_type) -> List:
    return list(map(to_type, inp))


def str_list_to_int_list(inp: List[str]) -> List[int]:
    return _str_list_to(inp, int)


def str_list_to_float_list(inp: List[str]) -> List[float]:
    return _str_list_to(inp, float)
