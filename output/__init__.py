#!/usr/bin/env python3
# created at Jul 23, 2017 4:09 PM by Nil-Zil

import sys

__all__ = ["eos", "read_file", "plot_check", "generate_test"]

pyversion = sys.executable
if 'python3' not in pyversion:
    print('Please use Python3 and above!')
