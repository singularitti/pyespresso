#!/usr/bin/env python3
# created at Jul 21, 2017 12:30 by Qi Zhang
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


class VCRelax:
    def __init__(self):
        self.input_dat = 'submit.dat'
        self.job_header = 'job_header'
        self.job_lines, self.schedule, self.modules, self.num_nodes, self.num_processors = self.read_job_header(
            self.job_header)
        self.pressures, self.v0, self.k0, self.k0p, self.vecs, self.vol_sc = self.read_input_params(
            self.input_dat)
        self.num_pressures = self.pressures.size
        self.error = 'CRASH'
        self.ps_cg, self.vs_cg, self.eos_opt_cg, self.eos_cov_cg = self._fit_crude_guess()

    @staticmethod
    def read_input_params(filename: str) -> tuple:
        """
        submit.dat has information on the calculation:
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
                'At least 7 are necessary for a good quality EOS fitting, please, modify your submit files.')
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
            for i in range(len(lines)):
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
            print('Therefore, I will consider' +
                  str(div) + 'processors per pressure.')
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
            print('Now I will generate the files for a crude guess scf calculation.')
            flow.create_files(uncalculated, (self.v0, self.k0,
                                             self.k0p), self.vecs, self.vol_sc, 'cg_', 'scf')
            flow.create_job(self.job_lines, uncalculated, self._distribute_task(), 'job-cg.sh', 'cg_', 'pw',
                            self.modules)

            job_sub = subprocess.Popen(
                ['sbatch', 'job-cg.sh'], stdout=subprocess.PIPE)
            output = job_sub.stdout.read().split()[-1]
            job_id = int(output)
            print('Job submitted. The job ID is' + str(job_id))
            print('Waiting calculation. This can take a while, you can grab a coffee!')

            # This part monitors the job execution.
            flag_job_cg = True
            while flag_job_cg:
                job_status = subprocess.Popen(['squeue', '-j', str(job_id)],
                                              stdout=subprocess.PIPE)  # request SLURM queue information
                outqueue = job_status.stdout.read().split()
                # read_file is a string that contains the job ID. If it is not in the queue information, the job is done.
                if output not in outqueue:
                    print(
                        "Crude guess is done! Continuing calculation... I hope your coffee was good!")
                    flag_job_cg = False
                else:  # If it is, sleeps a little bit and request queue information again.
                    time.sleep(2)
            flag_cg = False

    def read_crude_guess_output(self) -> tuple:
        """
        When crude guess step is done, we check if there are 'CRASH' files in the folder.
        And then we read the '.out' file from scf calculation, gain a list of P and a list of V.
        :return: (list, list)
        """
        num_error_files = 0
        ps = []
        vs = []

        for p in self.pressures:
            if self.error in os.listdir('cg_' + str(p)):
                # If 'CRASH' files exist, something went wrong. Stop the workflow.
                raise FileExistsError(
                    'Something went wrong, CRASH files found. Check your cg outputs.')
            else:
                output_file = 'cg_' + str(p) + '/' + 'cg_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as cg_out:
                        lines = cg_out.readlines()
                except FileNotFoundError:  # Check if read_file files exist.
                    print('The read_file file for P =' + str(p) +
                          'was not found. Removing it from list and continuing calculation.')
                    num_error_files += 1
                else:
                    for i in range(len(lines)):
                        if re.findall('P=', lines[i]):
                            ps.append(
                                float(lines[i].split('=')[-1].strip()) / 10)
                        if re.findall('volume', lines[i]):
                            vs.append(float(lines[i].split()[-2].strip()))

        if num_error_files > int(self.num_pressures / 2):
            # If more that a half of the calculations are not found, stop the workflow.
            raise ValueError(
                'More than a half of pressures were not calculated.')

        if set(ps) != set(self.pressures) or len(ps) != len(vs):
            # If the number of volumes are not equal to the number of pressures, something is wrong, stop the workflow.
            raise ValueError(
                'ERROR: Some pressures were not found in the files. Please, review your data.')

        return ps, vs

    def _fit_crude_guess(self) -> tuple:
        """
        This will use the submit V0, K0, and K0', as well as the list of P and V from read_crude_guess_output,
        to fit the Vinet equation of state.
        eos_opt are the returned V0, K0, and K0' for Vinet EOS.
        eos_cov are the estimated covariance of eos_opt.
        :return: (list, list, np.ndarray, np.ndarray)
        """
        ps, vs = self.read_crude_guess_output()
        eos_opt, eos_cov = curve_fit(
            flow.vinet, vs, ps, self.v0, self.k0, self.k0p)
        return ps, vs, eos_opt, eos_cov

    def write_crude_guess_result(self):
        """
        This subroutine creates a file named 'Initial_PxV.dat', and then writes the results fitted from the crude guess
        step.
        :return:
        """
        print('The initial volumes and pressures were written in the file Initial_PxV.dat')

        ps, vs, eos_opt, eos_cov = self._fit_crude_guess()
        std = np.sqrt(np.diag(eos_cov))

        with open('Initial_PxV.dat', 'w') as cg:
            cg.write('Initial PxV\n')
            cg.write('P (GPa)     V (au^3)\n')

            for i in range(0, len(ps)):
                cg.write("%7.2f   %7.2f\n" % (ps[i], vs[i]))

            print('EOS fitted, these are the standard deviations:')
            print('std V0 = %6.3f' % std[0])
            print('std K0 = %6.3f' % std[1])
            print('std Kp = %6.3f' % std[2])

            chi2 = 0
            for i in range(0, len(ps)):
                chi2 = chi2 + (ps[i] - flow.vinet(vs[i],
                                                  eos_opt[0], eos_opt[1], eos_opt[2])) ** 2
            chi = np.sqrt(chi2)
            print('chi = %7.4f' % chi)
            cg.write('Results for a Vinet EOS fitting:')
            cg.write('V0 = %6.4f    K0 = %4.2f    Kp = %4.2f' %
                     (eos_opt[0], eos_opt[1], eos_opt[2]))

    def vc_relax(self):
        """
        This subroutine creates 'vc_' folders, and then starts each vc-relaxation defined in each 'vc_' folders,
        and monitor it with a job_id.
        :return:
        """
        ps, vs, eos_opt, eos_cov = self._fit_crude_guess()

        vc_found = []

        for ps in self.pressures:
            if os.path.exists('vc_' + str(ps)):
                # Checks if cg_p folders exist.
                print('Folders vc for P = ' + str(ps) + 'found.')
                vc_found.append(ps)

        # Return the sorted, unique values in self.pressures that are not in cg_found.
        uncalculated: np.ndarray = np.setdiff1d(self.pressures, vc_found)
        # If flag_cg is True, continuously calculate.
        if uncalculated == np.array([]):
            print('All vc calculations exist. Proceeding to EOS fitting.')
            flag_vc = False
        else:
            print('There are remaining pressures to be calculated.')
            flag_vc = True

        # TODO: Here should be improved: If too many .out not found, do not rewrite .in, just use already existed files.
        while flag_vc:
            print('Now I will generate the files for a vc-relax calculation.')
            flow.create_files(uncalculated, eos_opt, self.vecs,
                              self.vol_sc, 'vc_', 'vc-relax')
            flow.create_job(self.job_lines, uncalculated, self._distribute_task(), 'job-vc.sh', 'vc_', 'pw',
                            self.modules)

            job_sub_vc = subprocess.Popen(
                ['sbatch', 'job-vc.sh'], stdout=subprocess.PIPE)
            output_vc = job_sub_vc.stdout.read().split()[-1]
            job_id_vc = int(output_vc)

            print(
                'Job for variable cell relaxation submitted. The job ID is' + str(job_id_vc))
            print(
                'Waiting calculation. This can take a long time, you can have lunch! Or take a nap! :)')

            # This part monitors the job execution.
            flag_job_vc = True
            while flag_job_vc:
                job_status_vc = subprocess.Popen(
                    ['squeue', '-j', str(job_id_vc)], stdout=subprocess.PIPE)
                outqueue_vc = job_status_vc.stdout.read().split()
                if output_vc not in outqueue_vc:
                    print("Congratulations, the variable cell relaxation is done!")
                    print("Fitting the final equation of state.")
                    flag_job_vc = False
                else:
                    time.sleep(2)
            flag_vc = False

    def read_vc_relax_output(self) -> tuple:
        """
        When vc-relax step is done, we check if there are 'CRASH' files in the folder.
        And then we read the '.out' file from vc-relax calculation, gain a list of final P and a list of V.
        :return: (list, list)
        """
        num_error_files = 0
        ps = []
        paux = []
        vs = []

        for p in self.pressures:
            if self.error in os.listdir('vc_' + str(p)):
                print(
                    'CRASH file found. Something went wrong with the vc-relax calculation, check your files.')
            else:
                output_file = 'vc_' + str(p) + '/' + 'vc_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as vc_out:
                        lines = vc_out.readlines()
                except FileNotFoundError:  # Check if read_file files exist.
                    print('The read_file file for P =' + str(p) +
                          'was not found. Removing it from list and continuing calculation.')  # This is a bug because you have less than 6 but it still works.
                    num_error_files += 1
                else:
                    for i in range(len(lines)):
                        if re.findall('P=', lines[i]):
                            paux.append(
                                float(lines[i].split('=')[-1].strip()) / 10)
                        if re.findall('volume', lines[i]):
                            vs.append(float(lines[i].split()[-2].strip()))
                    ps.append(paux[-1])

        return ps, vs

    def _fit_vc_relax(self):
        """
        This will use the submit V0, K0, and K0', as well as the list of P and V from read_vc_relax_output,
        to fit the Vinet equation of state.
        eos_opt are the returned V0, K0, and K0' for Vinet EOS.
        eos_cov are the estimated covariance of eos_opt.
        :return: (list, list, np.ndarray, np.ndarray)
        """
        ps, vs = self.read_vc_relax_output()
        eos_opt, eos_cov = curve_fit(flow.vinet, vs, ps, *self.eos_opt_cg)
        return ps, vs, eos_opt, eos_cov

    def write_vc_relax_result(self):
        ps, vs, eos_opt, eos_cov = self._fit_vc_relax()

        if len(ps) != len(vs):
            raise ValueError(
                'Some volumes are mossing, check your vc-relax outputs.')

        std = np.sqrt(np.diag(eos_cov))
        chi2 = 0
        for i in range(0, len(ps)):
            chi2 = chi2 + (ps[i] - flow.vinet(vs[i], *eos_opt)) ** 2

        with open('Out_file.dat', 'w') as vc:
            vc.write('P (GPa)    V (au^3)\n')
            for i in range(len(ps)):
                vc.write("%7.2f   %7.2f\n" % (ps[i], vs[i]))
                vc.write('Final Results for a Vinet EOS fitting:\n')
                vc.write('V0 = %6.4f    K0 = %4.2f    Kp = %4.2f' %
                         (eos_opt[0], eos_opt[1], eos_opt[2]))

        chi = np.sqrt(chi2)
        print('Ready, the vc-relax fitting is done. The standart deviations are:')
        print('V0 = %6.3f' % std[0])
        print('K0 = %6.3f' % std[1])
        print('Kp = %6.3f' % std[2])

        print('chi = %7.4f' % chi)
        print('The final EOS is written in the file Out_file.dat.')
        print('I hope you have a great day! See you next time.')
