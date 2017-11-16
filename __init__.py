#!/usr/bin/env python3
# created at Aug 26, 2017 2:56 PM by Qi Zhang

import sys

python_version = sys.executable
if 'python3' not in python_version:
    raise EnvironmentError('Please use Python executable higher than version 3.3!')

__author__ = {'Michel Lacerda': 'mld2189@columbia.edu', 'Qi Zhang': 'qz2280@columbia.edu',
              'Zhen Zhang': 'zz2427@columbia.edu'}
__copyright__ = 'Copyright 2017, Renata group'
__maintainer__ = 'Qi Zhang'
__credits__ = {'Renata M. M. Wentzcovitch': 'rmw2150@columbia.edu'}
__date__ = 'Nov 15, 2017'
