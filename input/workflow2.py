#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 21, 2017 12:30 by Nil-Zil

import os
import re

import numpy as np

print('Starting workflow calculations.')


class VCRelax(object):
    def __init__(self):
        self.inputdat = 'input.dat'
        self.jobheader = 'job_header'
        self.schedule, self.modules, self.nodes, self.processors = self.read_job_header(self.jobheader)
        self.pressures, self.v0, self.k0, self.k0p = self.read_input_params(self.inputdat)

    @staticmethod
    def read_input_params(filename) -> tuple:
        """
        input.dat has information on the calculation:
        The pressures to be calculated,
        the initial guess for V0, K0 and K0',
        and the lattice vectors, necessary to calculate the lattice parameter with the volume.
        :param filename: str
        :return: (np.ndarray, float, float, float)
        """
        with open(filename, 'r') as input_params:
            lines = input_params.readlines()

        pressures = np.array(lines[0].split(), dtype=float)
        num_pressures = pressures.size
        v0, k0, k0p = np.array(lines[1].split(), dtype=float)
        # lattice vectors a1, a2, a3
        a1 = np.array(lines[2].split(), dtype=float)
        a2 = np.array(lines[3].split(), dtype=float)
        a3 = np.array(lines[4].split(), dtype=float)
        vol_sc = np.fabs(np.dot(a3, np.cross(a1, a2)))  # volume of a simple cubic

        if num_pressures < 7:
            print('You have' + str(num_pressures) + 'pressure points!')
            print('At least 7 are necessary for a good quality EoS fitting, please, modify your input files.')
        else:
            print('You have' + str(num_pressures) + 'pressure points, more than 7! You are a smart user!')

        return pressures, v0, k0, k0p

    @staticmethod
    def read_job_header(filename) -> tuple:
        """
        job_header is a file with the parameters of the job to be submitted:
        wall time,
        number of nodes,
        number of processors,
        :param filename: str
        :return: tuple
        """
        with open(filename, 'r') as job_in:
            lines = job_in.readlines()
            for i in range(0, len(lines)):
                if re.findall('^Scheduler', lines[i]):
                    schedule = lines[i].split(':')[1].strip()
                if re.findall('^Necessary modules to be', lines[i]):
                    modules = lines[i].split(':')
                if re.findall('Number of nodes', lines[i]):
                    nodes = int(lines[i].split(':')[1].strip())
                if re.findall('Number of processor', lines[i]):
                    processors = int(lines[i].split(':')[1].strip())

        return schedule, modules, nodes, processors

    def _distribute_task(self):
        """
        Distribute calculations on different pressures evenly to all processors you have.
        :return: int
        """
        div = self.nodes * self.processors / self.pressures.size

        if (self.nodes * self.processors) % self.pressures.size == 0:
            print('Therefore, I will consider' + str(div) + 'procs per pressure.')
        else:
            raise ValueError('Number of processors not divided by number of pressures!')

        return int(div)

    def crude_guess(self):
        flag_cg = True
        num_cg = 0
        cg_found_list = []
        cont_cg = 0

        for j in self.pressures:
            if os.path.exists('cg_' + str(j)):
                print('Folders cg for P = ' + str(j) + 'found.')  # Checks if cg_j folders exist.
                cg_found_list.append(cont_cg)
                num_cg += 1
            cont_cg += 1

        if num_cg == self.pressures.size:
            print('All cg calculations exist. Proceeding to vc-relax calculation.')
            flag_cg = False
        elif num_cg == 0:
            print('Cg folders not found. Proceeding to crude guess calculation.')
            flag_cg = True
            press_cg = np.copy(self.pressures)
        else:
            print('Some remaining pressures to be calculated.')
            flag_cg = True
            press_cg = np.delete(self.pressures, cg_found_list)

        while flag_cg:
            functions.create_files(press_cg, eos_par, vecs, vol_sc, 'cg_', 'scf')
