#!/usr/bin/env python3
# created at Oct 20, 2017 6:13 PM by Qi Zhang

from pyque.meta.text import TextStream
from pyque.miscellaneous.strings import *

# ================================= These are some type aliases or type definitions. =================================
Namelist = TypeVar('Namelist')


def is_namelist(obj: object):
    if hasattr(obj, 'names') and hasattr(obj, 'default_values') and hasattr(obj, 'value_types'):
        return True
    else:
        return False


class SimpleParser(TextStream):
    def _match_one_string(self, pattern: str, *args):
        pass

    def _match_one_pattern(self, pattern: str, *args, **kwargs):
        """
        This method matches a pattern which exists only once in the file.

        :param pattern: a regular expression that you want to match
        :param args: a wrapper function which determines the returned type of value
        :return: Determined by the `wrapper`, the value you want to grep out from pyque.the file.
        """
        with open(self.infile, 'r') as f:
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
        This method reads each line simply from pyque.a file, and discards the '\n' character on the line end.

        :return: a list contains all the lines of the file
        """
        line_list = []
        with open(self.infile, 'r') as f:
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
        with open(self.infile, 'r') as f:
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
        This method reads 2 columns from pyque.a file. The are specified by spaces.

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
        with open(self.infile, 'r', encoding='utf-8') as f:
            for line in f:
                sp = line.split()
                keys.append(sp[col_index])
                del sp[col_index]  # Remove the indexing column
                values.append(func(sp))
        return dict(zip(keys, values))


class NamelistParserGeneric:
    def __init__(self, infile: str, namelist: object):
        """
        Match card between card title and the following '/' character.

        :param infile:
        :param namelist: a card has many names defined by Quantum ESPRESSO
        """
        self.infile = infile
        if is_namelist(namelist):
            self.namelist: Namelist = namelist
        else:
            raise TypeError('{0} is not a namelist!'.format(namelist))

    @staticmethod
    def _section_with_bounds(infile, start_pattern, end_pattern) -> Iterator[str]:
        """
        Search in file for the contents between 2 patterns. Referenced from
        [here](https://stackoverflow.com/questions/11156259/how-to-grep-lines-between-two-patterns-in-a-big-file-with-python).

        :param infile: file to be read
        :param start_pattern: the pattern labels where the content is going to start, the line contain this pattern is
            ignored
        :param end_pattern: the pattern labels where the content is to an end
        :return: an iterator that can read the file
        """
        section_flag = False
        for line in infile:
            if re.match(start_pattern, line, re.IGNORECASE):
                section_flag = True
                line = infile.readline()  # If the line is the `start_pattern` itself, we do not parse this line
            if line.startswith(end_pattern):
                section_flag = False
            if section_flag:
                yield line

    def read_namelist(self) -> Dict[str, str]:
        """
        A generic method to read a namelist.
        Note you cannot write more than one parameter in each line!

        :return: a dictionary that stores the inputted information of the intended card
        """
        namelist_names = set(self.namelist.names)
        filled_namelist: dict = {}
        start_pattern = '&' + self.namelist.__name__.upper()

        with open(self.infile, 'r') as f:
            generator: Iterator[str] = self._section_with_bounds(f, start_pattern, '/')  # '/' separates each namelist
            for line in generator:  # Read each line in the namelist until '/'
                s: str = line.strip()
                # Use '=' as the delimiter, split the stripped line into a key and a value.
                k, v = s.split('=', maxsplit=1)
                k: str = k.strip()
                v: str = v.strip().rstrip(',').strip()  # Ignore trailing comma of the line
                # Some keys have numbers as their labels, like 'celldm(i)', where $i = 1, \ldots, 6$. So we neet to
                # separate them.
                if '(' in k:
                    # Only take the part before the first '('
                    k_prefix = re.match("(\w+)\(?(\d*)\)?", k, flags=re.IGNORECASE).group(1)
                else:
                    k_prefix = k
                if k_prefix in namelist_names:
                    filled_namelist.update({k: v})
                else:
                    raise KeyError("'{0}' is not a valid name in '{1}' namelist!".format(k, self.namelist.__name__))

        return filled_namelist
