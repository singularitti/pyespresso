#!/usr/bin/env python3
# created at Oct 20, 2017 6:15 PM by Qi Zhang

import collections
import os

from read_file.read_basic import *


class PWscfOutputReader:
    @staticmethod
    def read_total_energy(inp: str):
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            match = re.findall(
                "!\s+total\s+energy\s+=\s+(-?\d+\.\d+)", f.read())
        return str_list_to_float_list(match)

    KMesh = NamedTuple(
        'KMesh', [('k_grid', List[float]), ('k_shift', List[float])])

    @staticmethod
    def read_k_mesh(inp: str) -> KMesh:
        """

        :param inp:
        :return:
        """
        with open(inp, 'r') as f:
            for line in f:
                if re.search('K_POINTS', line, re.IGNORECASE):
                    sp = f.readline().split()
                    k_grid = str_list_to_float_list(sp[0:3])
                    k_shift = str_list_to_float_list(sp[3:7])
                    k_mesh = collections.namedtuple(
                        'k_mesh', ['k_grid', 'k_shift'])
                    return k_mesh(k_grid, k_shift)
                else:
                    raise ValueError(
                        "'K_POINTS' not found in your submit! Please check!")

    @staticmethod
    def read_k_coordinates(in_file: str, out_file: str, coordinate_system: Optional[str] = 'crystal'):
        """
        This method can be used to read how many k-points are involved in a PWscf calculation from a file
        outputted by pw.x. Here regular expression is used. Different version of Quantum ESPRESSO may need different
        version of regular expression. If you find a bug, please contact the author at qz2280@columbia.edu.

        :param in_file: input file, where the data will be read from
        :param out_file: output file, where the data read
        :param coordinate_system: Can be 'Cartesian' at any case (upper, lower, etc.); or 'crystal' at any case.
        :return: None
        """
        regexp = "k\(\s+\d+\) = \((\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\), wk =\s+(\d+\.\d+)"
        cart_flag = "cart. coord. in units 2pi/alat"
        cryst_flag = "cryst. coord."

        if coordinate_system.lower() in {'cartesian', 'cart.', 'cart'}:
            system_type = cart_flag
        elif coordinate_system.lower() in {'crystal', 'cryst.', 'cryst'}:
            system_type = cryst_flag
        else:
            raise ValueError(
                "Unknown coordinate system type! It can be either 'Cartesian' or 'crystal'!")

        with open(in_file, 'r') as f, open(out_file, 'w') as g:
            flag = False  # If flag is true, read line and match pattern, if not true, no need to match pattern
            k_count = 0  # Count how many k-points have been read
            k_num = None  # How many k-points in total, given by Quantum ESPRESSO
            for line in f:
                if 'number of k points=' in line:
                    k_num = re.findall("number of k points=\s+(\d+)", line)[0]
                    g.write("{0}\n".format(k_num))
                if line.strip() == system_type:
                    flag = True
                if flag:
                    if k_count >= int(k_num):
                        flag = False

                    matches = re.finditer(regexp, line)
                    for match in matches:
                        # The first group is the k-point's coordinates 3-dimensional coordinates, and second
                        # is its weight in first Brillouin zone.
                        for group in match.groups():
                            g.write(group)
                            g.write('   ')  # A separation between first and second group
                        g.write("\n")

                    k_count += 1
        print('Reading done! File is stored at "{0}"'.format(os.path.abspath(out_file)))
