#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 21, 2017 12:30 by Nil-Zil
"""
This will do crude guess and vc-relax simulation automatically.
"""

import os
import re
import subprocess
import time

import numpy as np
from scipy.optimize import curve_fit

from . import flow

print('Starting workflow calculations.')


class VCRelax(object):
    def __init__(self):
        self.input_dat = 'input.dat'
        self.job_header = 'job_header'
        self.job_lines, self.schedule, self.modules, self.num_nodes, self.num_processors = self.read_job_header(
            self.job_header)
        self.pressures, self.v0, self.k0, self.k0p, self.vecs, self.vol_sc = self.read_input_params(
            self.input_dat)
        self.num_pressures = self.pressures.size
        self.error = 'CRASH'

    @staticmethod
    def read_input_params(filename: str) -> tuple:
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

        # Convert str to numpy.float64
        pressures = np.array(lines[0].split(), dtype=float)
        num_pressures = pressures.size
        if num_pressures < 7:
            print('You have' + str(num_pressures) + 'pressure points!')
            print(
                'At least 7 are necessary for a good quality EoS fitting, please, modify your input files.')
        else:
            print('You have' + str(num_pressures) +
                  'pressure points, more than 7! You are a smart user!')

        v0, k0, k0p = np.array(lines[1].split(), dtype=float)
        # lattice vectors a1, a2, a3
        a1 = np.array(lines[2].split(), dtype=float)
        a2 = np.array(lines[3].split(), dtype=float)
        a3 = np.array(lines[4].split(), dtype=float)
        vecs = np.stack((a1, a2, a3))
        # volume of a simple cubic
        vol_sc = np.fabs(np.dot(a3, np.cross(a1, a2)))

        return pressures, v0, k0, k0p, vecs, vol_sc

    @staticmethod
    def read_job_header(filename: str) -> tuple:
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
                    num_nodes = int(lines[i].split(':')[1].strip())
                if re.findall('Number of processor', lines[i]):
                    num_processors = int(lines[i].split(':')[1].strip())

        return lines, schedule, modules, num_nodes, num_processors

    def _distribute_task(self):
        """
        Distribute calculations on different pressures evenly to all processors you have.
        :return: int
        """
        print('You have' + str(self.num_pressures) + 'pressures and' + str(self.num_processors * self.num_nodes) +
              'processors.')
        div = self.num_nodes * self.num_processors / self.num_pressures

        if (self.num_nodes * self.num_processors) % self.num_pressures == 0:
            print('Therefore, I will consider' + str(div) + 'processors per pressure.')
        else:
            raise TypeError('Number of processors not divided by number of pressures!')

        return int(div)

    def crude_guess(self):
        """
        This subroutine creates 'cg_' folders, and then starts each crude guess defined in each 'cg_' folders,
        and monitor it with a job_id.
        :return:
        """
        cg_found = []

        for p in self.pressures:
            if os.path.exists('cg_' + str(p)):
                # Checks if cg_p folders exist.
                print('Folders cg for P = ' + str(p) + 'found.')
                cg_found.append(p)

        # Return the sorted, unique values in self.pressures that are not in cg_found.
        uncalculated: np.ndarray = np.setdiff1d(self.pressures, cg_found)
        # If flag_cg is True, continuously calculate.
        if uncalculated == np.array([]):
            print('All cg calculations exist. Proceeding to vc-relax calculation.')
            flag_cg = False
        else:
            print('There are remaining pressures to be calculated.')
            flag_cg = True

        while flag_cg:
            # Create jobs for those which have not been calculated.
            # If all have been calculated, skip this part and goes to vc-relax optimization.
            # The create_files accept an array, so it calculates the whole array, but I think maybe we can make it
            # accept one and loop against the array?
            flow.create_files(uncalculated, (self.v0, self.k0, self.k0p), self.vecs, self.vol_sc, 'cg_', 'scf')
            flow.create_job(self.job_lines, uncalculated, self._distribute_task(), 'job-cg.sh', 'cg_', 'pw',
                            self.modules)

            job_sub = subprocess.Popen(
                ['sbatch', 'job-cg.sh'], stdout=subprocess.PIPE)
            output = job_sub.stdout.read().split()[-1]
            job_id = int(output)
            print('Job submitted. The job ID is' + str(job_id))
            print('Waiting calculation. This can take a while, you can grab a coffee!')
            flag_job_cg = True

            # This part monitors the job execution.
            while flag_job_cg:
                job_status = subprocess.Popen(['squeue', '-j', str(job_id)],
                                              stdout=subprocess.PIPE)  # request SLURM queue information
                outqueue = job_status.stdout.read().split()
                # Output is a string that contains the job ID. If it is not in the queue information, the job is done.
                if output not in outqueue:
                    print("Crude guess is done! Continuing calculation... I hope your coffee was good!")
                    flag_job_cg = False
                else:  # If it is, sleeps a little bit and request queue information again.
                    time.sleep(2)
            flag_cg = False

    def check_crude_guess(self) -> tuple:
        """
        When crude guess step is done, we check if there are 'CRASH' files in the folder.
        :return: (list, list)
        """
        num_error_files = 0
        cont = 0
        p = []
        v = []

        for p in self.pressures:
            if self.error in os.listdir('cg_' + str(p)):
                # If 'CRASH' files exist, something went wrong. Stop the workflow.
                raise FileExistsError('Something went wrong, CRASH files found. Check your cg outputs.')
            else:
                output_file = 'cg_' + str(p) + '/' + 'cg_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as cg_out:
                        lines = cg_out.readlines()
                except FileNotFoundError:  # Check if output files exist.
                    print('The output file for P =' + str(p) +
                          'was not found. Removing it from list and continuing calculation.')
                    num_error_files += 1
                else:
                    for i in range(len(lines)):
                        if re.findall('P=', lines[i]):
                            p.append(float(lines[i].split('=')[-1].strip()) / 10)
                        if re.findall('volume', lines[i]):
                            v.append(float(lines[i].split()[-2].strip()))
                            # cont = cont + 1

        if num_error_files > int(self.num_pressures / 2):
            # If more that a half of the calculations are not found, stop the workflow.
            raise ValueError('More than a half of pressures were not calculated.')

        if set(p) != set(self.pressures) or len(p) != len(v):
            # If the number of volumes are not equal to the number of pressures, something is wrong, stop the workflow.
            raise ValueError('ERROR: Some pressures were not found in the files. Please, review your data.')

        return p, v

    def write_file(self):
        print('The initial volumes and pressures were written in the file Initial_PxV.dat')

        with open('Initial_PxV.dat', 'w') as pv:
            pv.write('Initial PxV\n')
            pv.write('P (GPa)     V (au^3)\n')
            p, v = self.check_crude_guess()
            for i in range(0, len(p)):
                pv.write("%7.2f   %7.2f\n" % (p[i], v[i]))
            eos_opt, eos_cov = curve_fit(flow.vinet, v, p, self.v0, self.k0, self.k0p)
            std = np.sqrt(eos_cov)
            print('EOS fitted, these are the standard deviations:')
            print('V0 = %6.3f' % std[0])
            print('K0 = %6.3f' % std[1])
            print('Kp = %6.3f' % std[2])

            chi2 = 0
            for i in range(0, len(p)):
                chi2 = chi2 + (p[i] - flow.vinet(v[i], eos_opt[0], eos_opt[1], eos_opt[2])) ** 2
            chi = np.sqrt(chi2)
            print('chi = %7.4f' % chi)
            pv.write('Results for a Vinet EoS fitting:')
            pv.write('V0 = %6.4f    K0 = %4.2f    Kp = %4.2f' % (eos_opt[0], eos_opt[1], eos_opt[2]))
