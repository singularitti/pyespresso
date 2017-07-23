#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 21, 2017 12:30 by Nil-Zil

print('Starting workflow calculations.')


class VCRelax(object):
    def __init__(self):
        self.input_dat = 'input.dat'
        self.job_header = 'job_header'
        self.job_lines, self.schedule, self.modules, self.nodes, self.processors = self.read_job_header(
            self.job_header)
        self.pressures, self.v0, self.k0, self.k0p, self.vecs, self.vol_sc = self.read_input_params(
            self.input_dat)

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
        vecs = np.stack((a1, a2, a3))
        # volume of a simple cubic
        vol_sc = np.fabs(np.dot(a3, np.cross(a1, a2)))

        if num_pressures < 7:
            print('You have' + str(num_pressures) + 'pressure points!')
            print(
                'At least 7 are necessary for a good quality EoS fitting, please, modify your input files.')
        else:
            print('You have' + str(num_pressures) +
                  'pressure points, more than 7! You are a smart user!')

        return pressures, v0, k0, k0p, vecs, vol_sc

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

        return lines, schedule, modules, nodes, processors

    def _distribute_task(self):
        """
        Distribute calculations on different pressures evenly to all processors you have.
        :return: int
        """
        div = self.nodes * self.processors / self.pressures.size

        if (self.nodes * self.processors) % self.pressures.size == 0:
            print('Therefore, I will consider' +
                  str(div) + 'procs per pressure.')
        else:
            raise ValueError(
                'Number of processors not divided by number of pressures!')

        return int(div)

    def crude_guess(self):
        """
        This subroutine creates 'cg_' folders, and then starts each crude guess defined in each 'cg_' folders,
        and monitor it with a job_id.
        :return:
        """
        num_cg = 0
        cg_found_list = []
        cont_cg = 0

        for p in self.pressures:
            if os.path.exists('cg_' + str(p)):
                # Checks if cg_p folders exist.
                print('Folders cg for P = ' + str(p) + 'found.')
                cg_found_list.append(cont_cg)
                num_cg += 1  # If found one, count one.
            cont_cg += 1  # No matter found or not, count one.

        # If flag_cg is True, continuously calculate.
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
            # Delete those which have been calculated.
            press_cg = np.delete(self.pressures, cg_found_list)

        while flag_cg:
            # Create jobs for those which have not been calculated.
            flow.create_files(press_cg, *(self.v0, self.k0, self.k0p), self.vecs, self.vol_sc, 'cg_', 'scf')
            flow.create_job(self.job_lines, press_cg, self._distribute_task(
            ), 'job-cg.sh', 'cg_', 'pw', self.modules)

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
                    print(
                        "Crude guess is done! Continuing calculation... I hope your coffee was good!")
                    flag_job_cg = False
                else:
                    time.sleep(2)

    def check_crude_guess(self):
        error = 'CRASH'
        num_error_files = 0
        cont = 0
        paux = []
        vaux = []

        for p in self.pressures:
            if error in os.listdir('cg_' + str(p)):
                # If CRASH files exist, something went wrong. Stop the workflow.
                raise FileExistsError(
                    'Something went wrong, CRASH files found. Check your cg outputs.')
            else:
                output_file = 'cg_' + str(p) + '/' + 'cg_' + str(p) + '.out'
                try:
                    with open(output_file, 'r') as qe_out:
                        qe_lines = qe_out.readlines()
                except FileNotFoundError:  # Check if output files exist.
                    print('The output file for P =' + str(p) +
                          'was not found. Removing it from list and continuing calculation.')
                    num_error_files += 1
                else:
                    for i in range(len(qe_lines)):
                        if re.findall('P=', qe_lines[i]):
                            paux.append(float(qe_lines[i].split(
                                '=')[-1].rstrip().lstrip()) / 10)
                        if re.findall('volume', qe_lines[i]):
                            vaux.append(
                                float(qe_lines[i].split()[-2].rstrip().lstrip()))
            cont = cont + 1

        if num_error_files > int(self.pressures.size / 2):
            # If more that a half of the calculations are not found, stop the workflow.
            raise ValueError(
                'More than a half of pressures were not calculated.')
