#!/usr/bin/env python3
# created at Jul 19, 2017 15:00 by Qi Zhang
"""
This module will only deal with elastic basics reading processes.
"""

from typing import *

import numpy as np

from pyque.parsers.simple import SimpleParser


class ElasticityOutputParser(SimpleParser):
    def read_elastic_tensor(self) -> Tuple[List[float], List[np.ndarray]]:
        """
        Read c_ij file, which looks like:
            Calculated tensor for each pressure:

            P = -10.0
             322.000    123.400    102.250      0.000      0.000      0.000
             123.400    322.000    102.250      0.000      0.000      0.000
             102.250    102.250    307.600      0.000      0.000      0.000
               0.000      0.000      0.000     72.175      0.000      0.000
               0.000      0.000      0.000      0.000     72.175      0.000
               0.000      0.000      0.000      0.000      0.000     99.300

            P = 0.0
             ...

        :return: A list of pressure and a list of 6 by 6 numpy arrays, each pressure corresponds to an array (tensor
            matrix).
        """
        pressures = []
        elastic_tensors = []
        with open(self.infile, 'r') as f:
            for line in f:
                if not line.strip():  # Ignore blank lines
                    continue
                if 'Calculated tensor for each pressure' in line:
                    f.readline()
                if 'P = ' in line:
                    pressures.append(float(line.split('=')[1]))
                    elastic_tensor = np.empty([6, 6], dtype=np.float64)  # Elastic tensor is a 6x6 matrix.
                    for i in range(6):
                        elastic_tensor[i] = list(map(float, f.readline().split()))
                    elastic_tensors.append(elastic_tensor)

        if len(pressures) == len(elastic_tensors):  # For each pressure there should be one elastic tensor
            return pressures, elastic_tensors
        else:
            raise ValueError(
                'There are {0} pressures, but {1} elastic tensors. Something went wrong, check your input file!'.format(
                    len(pressures), len(elastic_tensors)))

    def read_cij_vs_pressures(self) -> Dict[str, List[float]]:
        """
        Read first column as pressures, then the second to the last column as the compliance of the crystal at
        corresponding pressure.

        :return: a dictionary that contains pressure as one key, and a list of different elastic tensor's components as
            corresponding value.
        """
        return self._read_one_column_as_keys(0, lambda x: list(map(float, x)))
