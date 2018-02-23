#!/usr/bin/env python3
"""
This will do crude guess and vc-relax simulation automatically.
"""

import copy
import os
import pathlib
import re
import subprocess
import time
from typing import *
import schedule

import numpy as np

from pyque.core.qe_input import PWscfInput
from pyque.core.submitter import Submitter
from pyque.lexer.pwscf import PWscfOutputLexer
from pyque.tools.eos import VinetEoS


class VCRelaxSubmitter:
    def __init__(self, volumes, pwscf_inp, batch_inp):
        self.volumes = volumes
        self.pwscf_inp: PWscfInput = pwscf_inp
        self.batch_inp = batch_inp
        self.error = 'CRASH'

    @property
    def scaling_factors(self) -> List[float]:
        a1, a2, a3 = self.pwscf_inp.cell_parameters['value']
        v0 = np.fabs(np.dot(a3, np.cross(a1, a2)))
        return [(vol / v0) ** (1 / 3) for vol in self.volumes]

    def distribute_task(self) -> int:
        """
        Distribute calculations on different pressures evenly to all processors you have.

        :return: number of cores per pressure
        """
        scheduler = self.batch_inp
        tasks_number = len(self.volumes)
        cpus_per_task, remainder = divmod(scheduler.nodes_number * scheduler.cpus_per_node, tasks_number)
        if not remainder:  # If *remainder* is 0.
            print('Therefore, I will consider {0} cpus per volume.'.format(cpus_per_task))
        else:
            import warnings
            warnings.warn('Number of cpus is not divided by number of volumes!', stacklevel=2)
        return cpus_per_task

    def generate_crude_guess_folder(self):
        """
        This subroutine creates 'cg_' folders, and then starts each crude guess defined in each 'cg_' folders,
        and monitor it with a job_id.

        :return:
        """
        cg_found = set()

        for v in self.volumes:
            job_dir = 'cg_{0}'.format(v)
            if pathlib.Path(job_dir).is_dir():
                # Check if cg_ folders exist.
                print('Folders cg for V = {0} found.'.format(v))
                if pathlib.Path(job_dir + '/job.sh').exists() and pathlib.Path(job_dir + '/pwscf.in').exists():
                    cg_found.add(v)

        # Return the sorted, unique values in *self.volumes* that are not in cg_found.
        not_calculated: Set[float] = set(self.volumes) - cg_found

        if not_calculated == set():
            print('All cg calculations exist. Proceeding to vc-relax calculation.')
            return schedule.CancelJob
        else:
            print('There are remaining volumes {0} to be calculated.'.format(not_calculated))
            # Create jobs for those which have not been calculated.
            # If all have been calculated, skip this part and goes to vc-relax optimization.
            print('Now I will generate the folders for a crude guess scf calculation.')
            for i, v in enumerate(not_calculated):
                pathlib.Path('./cg_{0}'.format(v)).mkdir(parents=False, exist_ok=False)
                inp = copy.deepcopy(self.pwscf_inp)
                inp.system_namelist['celldm(1)'] *= self.scaling_factors[i]
                inp.to_text_file('./cg_{0}/pwscf.in'.format(v))
                self.batch_inp.directive.cpus_per_task = self.distribute_task()
                self.batch_inp.to_text_file('./cg_{0}/job.sh'.format(v))

    def submit_crude_guess(self):
        schedule.every(3).seconds.do(self.generate_crude_guess_folder)
        ids = []
        for v in self.volumes:
            path = 'cg_' + str(v)
            while not pathlib.Path(path + '/job.sh').exists() or not pathlib.Path(path + '/inp' + v + '.in').exists():
                # create job.sh or inp.in
                new_id = Submitter('job.sh').submit()
                if new_id:
                    ids.append(new_id)
                else:
                    time.sleep(5)
        if len(ids) == len(self.volumes):
            return ids
        else:
            return None

    def read_crude_guess_output(self) -> tuple:
        """
        When crude guess step is done, we check if there are 'CRASH' basics in the folder.
        And then we read the '.out' file from scf calculation, gain a list of P and a list of V.

        :return: (list, list)
        """
        error_files_num = 0
        ps = []
        vs = []

        for p in self.volumes:
            if self.error in os.listdir('cg_' + str(p)):
                # If 'CRASH' file exists, something went wrong. Stop the workflow.
                raise RuntimeError('Something went wrong, CRASH basics found. Check your cg outputs.')
            else:
                output_file = 'cg_' + str(p) + '/cg_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as cg_out:
                        lines = cg_out.read()
                except FileNotFoundError:  # Check if readers basics exist.
                    print('The file for P = {0} was not found. Remove it from list and continue calculation.'.format(p))
                    error_files_num += 1
                else:
                    lexer = PWscfOutputLexer(inp=lines)
                    vs.append(lexer.lex_cell_volume())
                    ps.append(lexer.lex_pressure())

        if set(ps) != set(self.volumes) or len(ps) != len(vs):
            # If the number of volumes are not equal to the number of pressures, something is wrong, stop the workflow.
            raise RuntimeError('Some pressures were not found in the basics. Please, review your data.')

        return ps, vs

    def write_crude_guess_result(self):
        """
        This subroutine creates a file named 'Initial_PxV.dat', and then writes the results fitted from the crude guess
        step.

        :return:
        """
        print('The initial volumes and pressures were written in the file Initial_P_vs_V.')

        ps, vs = self.read_crude_guess_output()
        eos_opt, eos_cov = self.vinet.fit_p_of_v(vs, ps)
        std = np.sqrt(np.diag(eos_cov))

        with open('Initial_P_vs_V', 'w') as f:
            f.write('P (GPa)     V (au^3)\n')

            for i in range(0, len(ps)):
                f.write("{:7.2f}   {:7.2f}\n".format(ps[i], vs[i]))

            print('EOS fitted, these are the standard deviations:')
            print('std V0 = {:6.3f}'.format(std[0]))
            print('std K0 = {:6.3f}'.format(std[1]))
            print('std K0p = {:6.3f}'.format(std[2]))

            chi2 = 0
            for i in range(0, len(ps)):
                chi2 = chi2 + (ps[i] - VinetEoS(eos_opt[0], eos_opt[1], eos_opt[2]).p_of_v(vs[i])) ** 2
            chi = np.sqrt(chi2)
            print('chi = {:7.4f}'.format(chi))
            f.write('Results for a VinetEoS EOS fitting:')
            f.write('V0 = {:6.4f}    K0 = {:4.2f}    K0p = {:4.2f}'.format(eos_opt[0], eos_opt[1], eos_opt[2]))

    def vc_relax(self):
        """
        This subroutine creates 'vc_' folders, and then starts each vc-relaxation defined in each 'vc_' folders,
        and monitor it with a job_id.

        :return:
        """
        ps, vs, eos_opt, eos_cov = self._fit_crude_guess()

        vc_found = []

        for ps in self.volumes:
            if os.path.exists('vc_' + str(ps)):
                # Checks if cg_p folders exist.
                print('Folders vc for P = ' + str(ps) + 'found.')
                vc_found.append(ps)

        # Return the sorted, unique values in self.pressures that are not in cg_found.
        uncalculated: np.ndarray = np.setdiff1d(self.volumes, vc_found)
        # If flag_cg is True, continuously calculate.
        if uncalculated == np.array([]):
            print('All vc calculations exist. Proceeding to EOS fitting.')
            flag_vc = False
        else:
            print('There are remaining pressures to be calculated.')
            flag_vc = True

        # TODO: Here should be improved: If too many .out not found, do not rewrite .in, just use already existed basics.
        while flag_vc:
            print('Now I will generate the basics for a vc-relax calculation.')
            flow.create_files(uncalculated, eos_opt, self.vecs,
                              self.vol_sc, 'vc_', 'vc-relax')
            flow.create_job(self.job_lines, uncalculated, self.distribute_task(), 'job-vc.sh', 'vc_', 'pw',
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
        When vc-relax step is done, we check if there are 'CRASH' basics in the folder.
        And then we read the '.out' file from vc-relax calculation, gain a list of final P and a list of V.

        :return: (list, list)
        """
        num_error_files = 0
        ps = []
        paux = []
        vs = []

        for p in self.volumes:
            if self.error in os.listdir('vc_' + str(p)):
                print(
                    'CRASH file found. Something went wrong with the vc-relax calculation, check your basics.')
            else:
                output_file = 'vc_' + str(p) + '/' + 'vc_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as vc_out:
                        lines = vc_out.readlines()
                except FileNotFoundError:  # Check if readers basics exist.
                    print('The readers file for P =' + str(p) +
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

    def write_vc_relax_result(self):
        ps, vs = self.read_vc_relax_output()
        eos_opt, eos_cov = self.vinet.fit_p_of_v(vs, ps)

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
                vc.write('Final Results for a VinetEoS EOS fitting:\n')
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
