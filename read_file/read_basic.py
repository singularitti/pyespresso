#!/usr/bin/env python3
# created at Oct 20, 2017 6:13 PM by Qi Zhang

import re
from typing import *


class SimpleRead:
    @staticmethod
    def read_each_line(file: str) -> List[str]:
        """
        This method reads each line simply from a file, and discard the '\n' character.

        :param file: A single file that is to be read.
        :return: A list contains all the lines of the file.
        """
        line_list = []
        with open(file, 'r') as f:
            for line in f:
                line_list.append(re.split('\n', line)[0])
        return line_list

    @staticmethod
    def read_two_columns(file: str) -> Tuple[List[str], List[str]]:
        """
        This method reads 2 columns from a file.

        :param file: A single file that is to be read.
        :return: two columns of the file
        """
        col1_list = []
        col2_list = []
        with open(file, 'r') as f:
            for line in f:
                if not line.split():
                    f.readline()
                else:
                    sp = line.split()
                    col1_list.append(sp[0])
                    col2_list.append(sp[1])
        return col1_list, col2_list

    @staticmethod
    def read_one_column_as_keys(file: str, col_index: int, wrapper: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param file: A single file that is to be read.
        :param col_index: the index of the column that you want to make it as keys.
        :param wrapper: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        key_list = []
        value_list = []
        # Add utf-8 support because we may use special characters.
        with open(file, 'r', encoding='utf-8') as f:
            for line in f:
                sp = line.split()
                key_list.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                value_list.append(wrapper(sp))
        return dict(zip(key_list, value_list))

    def _read_reciprocal_points(self, file: str) -> Dict[str, List[float]]:
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

        :param file: file you want to specify your k-points
        :return: a dictionary
        """
        return self.read_one_column_as_keys(file, 0, lambda x: list(map(float, x)))


def _str_list_to(inp: List[str], to_type) -> List:
    return list(map(to_type, inp))


def str_list_to_int_list(inp: List[str]) -> List[int]:
    return _str_list_to(inp, int)


def str_list_to_float_list(inp: List[str]) -> List[float]:
    return _str_list_to(inp, float)
