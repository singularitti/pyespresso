#! /usr/bin/env python3



working_dir = os.getcwd()
dirs_vc = os.listdir(working_dir)
erro = 'CRASH'

signal_exist = True
for i in dirs_vc:
	if i.startswith('vc_'):
#		print('Folders with vc-relax calculations found. Performing initial scf calculations.')
		signal_exist = False
		break
if signal_exist:
	print('Folders with vc-relax calculations were not found. You need to optimize the structure at each pressure before the phonon calculation. Exiting.')
	quit()

with open('input.dat') as input_params:
	lines_input = input_params.readlines()

press = np.array(lines_input[0].split(), dtype=float)
npress = len(press)
# Reads numnber of atoms and number of atomic types from input
with open('qe_input_data','r') as in_data:
	lines = in_data.readlines()
	for i in range(0, len(lines)):
		if lines[i].startswith('number of atoms'):
			nat = int(lines[i].split(':')[1].rstrip().lstrip())
		if lines[i].startswith('number of atomic'):
			ntyp = lines[i].split(':')[1].rstrip().lstrip()


with open('job_header','r') as job_in:
	job_lines = job_in.readlines()
	for i in range(0,len(job_lines)):
		if re.findall('^Scheduler', job_lines[i]):
			schedul = job_lines[i].split(':')[1].rstrip().lstrip()
		if re.findall('^Necessary modules to be', job_lines[i]):
			modul = job_lines[i].split(':')
		if re.findall('Number of nodes', job_lines[i]):
			nodes = int(job_lines[i].split(':')[1].rstrip().lstrip())
		if re.findall('Number of processor', job_lines[i]):
			procs = int(job_lines[i].split(':')[1].rstrip().lstrip())

div = (nodes*procs)/npress

if (nodes*procs) % npress == 0:
	print('I will consider %d procs per pressure.' % div)
else:
	print('Number of processors not divided by number of pressures. Exiting')
	quit()

#print('Generating structure files')

signal_scf = True

p_not_found_list = []
number_err_files = 0
cont = 0
atom_pos = []
for i in press:
	del atom_pos[:]
	if erro in os.listdir('vc_' + str(i)):
		print('CRASH files found. Something went wrong with the vc-relax calculation for P = %s. Check your vc_ folders.' % str(i))
		quit()
	else:
		file_path = 'vc_' + str(i) + '/' + 'vc_' + str(i) + '.out'
		try:
			with open(file_path) as vc_out:
				vc_lines = vc_out.readlines()
		except FileNotFoundError: # Check if vc output file exists
			print('The output file for P = %d was not found. Removing it from pressure list and continuing calculation.' % i)
			p_not_found_list = p_not_found_list + [cont]
			number_err_files = number_err_files + 1 # counts the number of files that wasn't found
			pass
		else:
			for j in range(0,len(vc_lines)):
				if vc_lines[j].startswith('Begin final'):
					alat = vc_lines[j-1].split()[3].rstrip().lstrip()
					a1 = np.array(vc_lines[j+5].split(), dtype=float)
					a2 = np.array(vc_lines[j+6].split(), dtype=float)
					a3 = np.array(vc_lines[j+7].split(), dtype=float)
					for k in range(0,nat):
						atom_pos = atom_pos + [vc_lines[j + 10 + k]]
	vectors = np.array([a1,a2,a3])
	dir_name = 'scf_' + str(i)
	phon_name = dir_name + '-ph.in'
	if os.path.exists(dir_name):
		print('scf folder for P = %d exists. Skiping qe input generation and creating phonon files.' % i)
	else:
		print('scf folder for P = %d does not exist. Creating input files for pw and ph.' % i)
		functions.create_scf_files(vectors, alat, dir_name, int(i), atom_pos)
	if os.path.exists(dir_name + '/' + dir_name + '.out'):
		print('Scf out files exist. Skiping scf calculation')
		signal_scf = False
	functions.create_ph_input(phon_name)
	shutil.move(phon_name, dir_name)


if signal_scf:
	signal_job_scf = True
	functions.create_job(job_lines, press, div, 'job-scf.sh', 'scf_', 'pw', modul)
	job_sub_scf = subprocess.Popen(['sbatch', 'job-scf.sh'], stdout=subprocess.PIPE)
	output_scf = job_sub_scf.stdout.read().split()[-1]
	job_id_scf = int(output_scf)
	
	print('Job submitted. The Job ID is %d' % job_id_scf)
	print('Waiting calculation. These are usually fast.')
	
	while signal_job_scf:
		job_status_scf = subprocess.Popen(['squeue', '-j', str(job_id_scf)], stdout=subprocess.PIPE)
		outqueue_scf = job_status_scf.stdout.read().split()
		if output_scf not in outqueue_scf:
			print("The scf calculation is done! Submitting phonon calculation.")
			signal_job_scf = False
		else:
			time.sleep(2)
	

for i in press:
	if erro in os.listdir('scf_' + str(i)):
		print('CRASH files found. Check your scf calculations')
		quit()

functions.create_job(job_lines, press, div, 'job-ph.sh', 'scf_', 'ph', modul)
job_sub_ph = subprocess.Popen(['sbatch', 'job-ph.sh'], stdout=subprocess.PIPE)
output_ph = job_sub_ph.stdout.read().split()[-1]
job_id_ph = int(output_ph)

print('Phonon Job submitted. The Job ID is %d' % job_id_ph)
print('These calculations take a long time. I will exit now, continue monitoring your job!')
print('Have a nice day!')
