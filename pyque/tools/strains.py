#!/usr/bin/env python3
"""
:mod:`strain` -- title
========================================

.. module strain
   :platform: Unix, Windows, Mac, Linux
   :synopsis: doc
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""


class EulerianStrain:
    def __init__(self, v0: float, v: float):
        self.v0 = v0
        self.v = v

    @property
    def isotropic_compression(self) -> float:
        return 1 / 2 * ((self.v0 / self.v) ** (2 / 3) - 1)

    def is_compression(self) -> bool:
        if self.isotropic_compression > 0:
            return True
        else:
            return False


class LagrangianStrain:
    def __init__(self, v0: float, v: float):
        self.v0 = v0
        self.v = v

    @property
    def isotropic_compression(self) -> float:
        return 1 / 2 * ((self.v / self.v0) ** (2 / 3) - 1)

    def is_compression(self) -> bool:
        if self.isotropic_compression < 0:
            return True
        else:
            return False
