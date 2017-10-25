#!/usr/bin/env python3
# created at Oct 24, 2017 6:17 PM by Qi Zhang
"""
Referenced from https://d32ogoqmya1dw8.cloudfront.net/images/NAGTWorkshops/mineralogy/mineral_physics/table_9.v13.png.
"""

HEXAGONAL = {'indices': [(0, 0), (2, 2), (3, 3), (0, 1), (0, 2)],
             'names': ['c11', 'c33', 'c44', 'c12', 'c31']}

ISOTROPIC = [(0, 0), (0, 1)]

CUBIC = [(0, 0), (3, 3), (0, 1)]

ORTHORHOMBIC = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (0, 2), (1, 2)]

crystal_classes = {'hex': HEXAGONAL, 'iso': ISOTROPIC, 'cubic': CUBIC, 'orth': ORTHORHOMBIC}
