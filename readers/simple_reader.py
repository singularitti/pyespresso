#!/usr/bin/env python3
# created at Oct 20, 2017 6:13 PM by Qi Zhang

from miscellaneous.string import *


class SingleFileReader:
    def __init__(self, in_file: str):
        """
        In our implementation, for each file, you need to generate a `SimpleRead` class to read it.

        :param in_file: the exact one input file for this class
        """
        self.in_file = in_file

    def _match_only_once(self, pattern: str, wrapper: Callable):
        """
        This method matches a pattern which exists only once in the file.

        :param pattern: a regular expression that you want to match
        :param wrapper: a wrapper function which determines the returned type of value
        :return: Determined by the `wrapper`, the value you want to grep out from the file.
        """
        with open(self.in_file, 'r') as f:
            match = re.findall(pattern, f.read())
        return wrapper(match[0])

    def read_line_by_line(self) -> List[str]:
        """
        This method reads each line simply from a file, and discards the '\n' character on the line end.

        :return: a list contains all the lines of the file
        """
        line_list = []
        with open(self.in_file, 'r') as f:
            for line in f:
                line_list.append(re.split('\n', line)[0])
        return line_list

    def _read_n_columns(self, n: int) -> List[List[str]]:
        """
        Read any number of columns (if exists) separated by spaces.

        :param n: an integer specifies the number of columns you want to read
        :return: a list of `n` columns that contain the contents of each column
        """
        n_columns = [[] for i in range(n)]
        with open(self.in_file, 'r') as f:
            for line in f:
                if not line.split():  # If line is ''
                    f.readline()
                else:
                    sp = line.split()
                    try:
                        for i, column in enumerate(n_columns):
                            column.append(sp[i])
                    except IndexError:
                        print('You want more columns than you have! Check your file at line {0}!'.format(line))
        return n_columns

    def read_two_columns(self) -> List[List[str]]:
        """
        This method reads 2 columns from a file. The are specified by spaces.

        :return: two columns of the file
        """
        return self._read_n_columns(2)

    def _read_one_column_as_keys(self, col_index: int, func: Callable[[List[str]], Any]) -> Dict[str, Any]:
        """
        This method read the one of the columns of a file as keys,
        the combination of rest columns are values to corresponding keys.

        :param col_index: the index of the column that you want to make it as keys.
        :param func: A function that can process the values to the form that you want.
        :return: A dictionary that could contain anything as its values, but with strings as its keys.
        """
        keys = []
        values = []
        # Add utf-8 support because we may use special characters.
        with open(self.in_file, 'r', encoding='utf-8') as f:
            for line in f:
                sp = line.split()
                keys.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                values.append(func(sp))
        return dict(zip(keys, values))

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
        return self._read_one_column_as_keys(0, lambda x: list(map(float, x)))
