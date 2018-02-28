#!/usr/bin/env python3
# created at Dec 2, 2017 5:23 PM by Qi Zhang
"""
checkvcfin.sh
"""


def stream_has_colours(stream):
    """
    True if stream supports colours. Python cookbook, #475186
    """
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        return False  # auto color only on TTYs
    try:
        import curses
        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except:
        return False  # guess false in case of error


class StringColorizer(object):
    """
    From [here](https://raw.githubusercontent.com/materialsproject/pymatgen/master/pymatgen/util/string.py).
    """
    colours = {"default": "",
               "blue": "\x1b[01;34m",
               "cyan": "\x1b[01;36m",
               "green": "\x1b[01;32m",
               "red": "\x1b[01;31m",
               # lighting colours.
               # "lred":    "\x1b[01;05;37;41m"
               }

    def __init__(self, stream):
        self.has_colours = stream_has_colours(stream)

    def __call__(self, string, colour):
        if self.has_colours:
            code = self.colours.get(colour.lower(), "")
            if code:
                return code + string + "\x1b[00m"
            else:
                return string
        else:
            return string
