#!/usr/bin/env python3
# created at Jul 21, 2017 12:29 by Nil-Zil

import re


# TODO: Separate by :, [], ! comment

class SCFGenerator:
    def __init__(self):
        pass

    def readfile(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                kv = lines[i].split('=', maxsplit=1)
                k = kv[0]
                v = kv[1:]
                if '!' in v:
                    v = v.split('!')[0]  # ignore things after the first '!'
