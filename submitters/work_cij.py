#! /usr/bin/env python3

erro = 'CRASH'

with open('input_cij') as cij_in:
    lines_cij = cij_in.readlines()
    for i in range(0, len(lines_cij)):
        if lines_cij[i].startswith('Crystal Class'):
            cclass = lines_cij[i].split(':')[1].rstrip().lstrip()
        if lines_cij[i].startswith('Strain'):
            dist = []
            dist = lines_cij[i].split()
        if lines_cij[i].startswith('Submit'):
            job_will_be_submited = lines_cij[i].split(':')[1].rstrip().lstrip()

if job_will_be_submited == 'y' or job_will_be_submited == 'Y':
    submit_job = True
else:
    submit_job = False

nstrain = func_cij.define_cij(cclass)
number_dists = len(dist) - 1
print('You need %d strains.' % nstrain)

with open('submitters.dat') as input_params:
    lines_input = input_params.readlines()

press = np.array(lines_input[0].split(), dtype=float)
npress = len(press)

print('Checking previously Cij calculation.')
number_cij = 0
for i in press:
    if os.path.exists('scf_' + str(i) + '/' + 'cij'):
        print('Cij folder for P = %s exists.' % str(i))
        number_cij = number_cij + 1

if number_cij == npress:
    print('All pressures have their Cij folder. Proceeding to the tensor calculation.')
    calculate_tensor = True
elif number_cij == 0:
    print('Cij folders do not exist. Preparing files for their calculation.')
    calculate_tensor = False
else:
    print('Not all pressures have their Cij folder. Preparing files for the missing ones')


# Reads numnber of atoms and number of atomic types from submitters
with open('qe_input_data', 'r') as in_data:
    lines = in_data.readlines()
    for i in range(0, len(lines)):
        if lines[i].startswith('number of atoms'):
            nat = int(lines[i].split(':')[1].rstrip().lstrip())
        if lines[i].startswith('number of atomic'):
            ntyp = lines[i].split(':')[1].rstrip().lstrip()

with open('job_header', 'r') as job_in:
    job_lines = job_in.readlines()
    for i in range(0, len(job_lines)):
        if re.findall('^Scheduler', job_lines[i]):
            schedul = job_lines[i].split(':')[1].rstrip().lstrip()
        if re.findall('^Necessary modules to be', job_lines[i]):
            modul = job_lines[i].split(':')
        if re.findall('Number of nodes', job_lines[i]):
            nodes = int(job_lines[i].split(':')[1].rstrip().lstrip())
        if re.findall('Number of processor', job_lines[i]):
            procs = int(job_lines[i].split(':')[1].rstrip().lstrip())

if not calculate_tensor:

    total_number_of_calculations = npress * nstrain * number_dists
    print(total_number_of_calculations)
    div = (nodes * procs) / total_number_of_calculations
    if (nodes * procs) % total_number_of_calculations == 0:
        div = (nodes * procs) / total_number_of_calculations
        print('I will consider %d procs per pressure.' % div)
    elif (nodes * procs) < total_number_of_calculations:
        print('You have more calculations than processors, considering breaking you job in more than one.')
        print('The job file will be created, but it will not be submitted')
        submit_job = False
    else:
        print('Number of processors not divided by number of pressures. Exiting')
        print('The job file will be created, but it will not be submitted')
        submit_job = False

    signal_scf = True
    p_not_found_list = []
    number_err_files = 0
    cont = 0
    atom_pos = []

    for i in press:
        del atom_pos[:]
        dir_name = 'scf_' + str(i)
        if erro in os.listdir('vc_' + str(i)):
            print('CRASH files found. Something went wrong with the vc-relax calculation for P = %s. Check your vc_ folders.' % str(i))
            quit()
        else:
            file_path = 'vc_' + str(i) + '/' + 'vc_' + str(i) + '.out'
            try:
                with open(file_path) as vc_out:
                    vc_lines = vc_out.readlines()
            except FileNotFoundError:  # Check if vc output file exists
                print(
                    'The readers file for P = %d was not found. Removing it from pressure list and continuing calculation.' % i)
                p_not_found_list = p_not_found_list + [cont]
                # counts the number of files that wasn't found
                number_err_files = number_err_files + 1
                pass
            else:
                for j in range(0, len(vc_lines)):
                    if vc_lines[j].startswith('Begin final'):
                        alat = vc_lines[j - 1].split()[3].rstrip().lstrip()
                        a1 = np.array(vc_lines[j + 5].split(), dtype=float)
                        a2 = np.array(vc_lines[j + 6].split(), dtype=float)
                        a3 = np.array(vc_lines[j + 7].split(), dtype=float)
                        for k in range(0, nat):
                            atom_pos = atom_pos + [vc_lines[j + 10 + k]]
        vectors = np.array([a1, a2, a3])
        if os.path.exists(dir_name):
            print(
                'Scf folder for P = %d exists. Skiping qe submitters generation and creating Cij files.' % i)
            signal_scf = False
        else:
            print(
                'Scf folder for P = %d does not exist. Creating submitters files for pw and Cij.' % i)
            functions.create_scf_files(
                vectors, alat, dir_name, int(i), atom_pos)
        cij_dir_name = 'cij'
        if os.path.exists(dir_name + '/' + dir_name + '.out'):
            print('Scf out files exist. Skiping scf calculation')
            signal_scf = False
        if os.path.exists(dir_name + '/' + cij_dir_name):
            print('Cij folder for pressure P = %s exists. ' % str(i))
        else:
            functions.create_cij_files(
                vectors, alat, dir_name, cij_dir_name, i, atom_pos, dist, nstrain, cclass)

    if signal_scf:
        div_scf = int((nodes * procs) / npress)
        if (nodes * procs) % npress == 0:
            print('Considering %d procs per pressure.' % div_scf)
            signal_job_scf = True
        else:
            print('Number of processors not divided by number of pressures. Job file will be created, but will not be submitted.')
            signal_job_scf = False
        functions.create_job(job_lines, press, div_scf,
                             'job-scf.sh', 'scf_', 'pw', modul)
        job_sub_scf = subprocess.Popen(
            ['sbatch', 'job-scf.sh'], stdout=subprocess.PIPE)
        output_scf = job_sub_scf.stdout.read().split()[-1]
        job_id_scf = int(output_scf)

        if signal_job_scf:
            print('Job submitted. The Job ID is %d' % job_id_scf)
            print('Waiting calculation. These are usually fast.')
        else:
            print('Starting Cij calculation without performing scf ones.')

        while signal_job_scf:
            job_status_scf = subprocess.Popen(
                ['squeue', '-j', str(job_id_scf)], stdout=subprocess.PIPE)
            outqueue_scf = job_status_scf.stdout.read().split()
            if output_scf not in outqueue_scf:
                print("The scf calculation is done! Preparing Cij calculation.")
                signal_job_scf = False
            else:
                time.sleep(2)

    functions.create_job_cij(
        job_lines, press, div, 'job-cij.sh', 'scf_', cclass, nstrain, dist, modul)

    if submit_job:
        signal_job_cij = True
        job_sub_cij = subprocess.Popen(
            ['sbatch', 'job-cij.sh'], stdout=subprocess.PIPE)
        output_cij = job_sub_cij.stdout.read().split()[-1]
        job_id_cij = int(output_cij)

        print('Job submitted. The Job ID is %d' % job_id_cij)
        print('After the Job is done, run this script again to calculate the Cij tensor for all pressures.')

    #	while signal_job_cij:
    #		job_status_cij = subprocess.Popen(['squeue', '-j', str(job_id_cij)], stdout=subprocess.PIPE)
    #		outqueue_cij = job_status_cij.stdout.read().split()
    #		if output_cij not in outqueue_cij:
    #			print("The scf calculation is done! Submitting phonon calculation.")
    #			signal_job_cij = False
    #		else:
    #			time.sleep(2)
cij_file_name = 'Cij_file.dat'

fp = open('CijxP.dat', 'w')

if calculate_tensor:
    cij_file = open(cij_file_name, 'w')
    cij_file.write("Calculated tensor for each pressure:\n")
    Cij = np.zeros((6, 6))
    for i in press:
        cij_file.write('\nP = %s\n' % str(i))
        dir_name = 'scf_' + str(i)
        dir_cij = dir_name + '/' + 'cij'
        for j in range(1, nstrain + 1):
            functions.calculate_cij_tensor(Cij, dir_cij, cclass, dist, j)
        functions.print_cij(cij_file, Cij)
        fp.write('%6.2f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n' %
                 (i, Cij[0, 0], Cij[2, 2], Cij[3, 3], Cij[0, 1], Cij[0, 2]))
    cij_file.close()
