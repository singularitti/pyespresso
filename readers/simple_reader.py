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

    def _match_one_string(self, pattern: str, *args):
        pass

    def _match_one_pattern(self, pattern: str, *args, **kwargs):
        """
        This method matches a pattern which exists only once in the file.

        :param pattern: a regular expression that you want to match
        :param args: a wrapper function which determines the returned type of value
        :return: Determined by the `wrapper`, the value you want to grep out from the file.
        """
        with open(self.in_file, 'r') as f:
            s = f.read()
            match: Optional[List[str]] = re.findall(pattern, s,
                                                    **kwargs)  # `match` is either an empty list or a list of strings.
            if match:
                if len(args) == 0:  # If no wrapper argument is given, return directly the matched string
                    return match
                elif len(args) == 1:  # If wrapper argument is given, i.e., not empty, then apply wrapper to the match
                    wrapper, = args
                    return [wrapper(m) for m in match]
                else:
                    raise TypeError('Multiple wrappers are given! Only one should be given!')
            else:  # If no match is found
                print('Pattern {0} not found in string {1}!'.format(pattern, s))
                return None

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
        n_columns = [[] for _ in range(n)]
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


class MultipleFilesReader:
    def __init__(self, files: List[str]):
        # Construct a dictionary with file names as keys
        self.files = {file: SingleFileReader(file) for file in files}
