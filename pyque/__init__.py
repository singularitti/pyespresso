#!/usr/bin/env python3

import sys

PY3 = sys.version_info

if PY3 < (3, 5):
    raise EnvironmentError('Please use Python version higher than 3.5!')

__author__ = {'Michel Lacerda': 'mld2189@columbia.edu', 'Qi Zhang': 'qz2280@columbia.edu'}
__copyright__ = 'Copyright (c) 2017, Renata group'
__credits__ = {'Renata M. M. Wentzcovitch': 'rmw2150@columbia.edu'}
__date__ = 'Nov 15, 2017'
__maintainer__ = 'Qi Zhang'
