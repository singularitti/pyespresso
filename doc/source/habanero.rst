Habanero tutorials
==================

Getting Started
---------------

See `Habanero "Getting Started" page <https://confluence.columbia.edu/confluence/display/rcs/Habanero+-+Getting+Started/>`_
for references.

This is a list of the available softwares on Habanero:

+---------+-----+------------------+------------------+--------------+
| Name    | Ver | Location /       | RPM / files      | Category     |
|         | sio | Module           |                  |              |
|         | n   |                  |                  |              |
+=========+=====+==================+==================+==============+
| anacond | pyt | module load      | /rigel/opt/anaco | Python for   |
| a-4.2.0 | hon | anaconda/2-4.2.0 | nda2-4.2.0       | Scientific   |
|         | 2.7 |                  |                  | Computing    |
|         | .12 |                  |                  |              |
+---------+-----+------------------+------------------+--------------+
| anacond | pyt | module load      | /rigel/opt/anaco | Python for   |
| a-4.2.0 | hon | anaconda/3-4.2.0 | nda3-4.2.0       | Scientific   |
|         | 3.5 |                  |                  | Computing    |
|         | .2  |                  |                  |              |
+---------+-----+------------------+------------------+--------------+
| beagle  | 2.1 | module load      | /rigel/opt/beagl | Genetic      |
|         | .2  | beagle/2.1.2     | e-2.1.2          | Analysis     |
+---------+-----+------------------+------------------+--------------+
| beast   | 1.8 | module load      | /rigel/opt/BEAST | Genetic      |
|         | .4  | beast/1.8.4      | v1.8.4           | Analysis     |
+---------+-----+------------------+------------------+--------------+
| cuda    | 8.0 | module load      | /cm/shared/apps/ | GPU          |
|         | .44 | cuda80/toolkit/8 | cuda80/toolkit/8 | Computing    |
|         |     | .0.44            | .0.44            |              |
+---------+-----+------------------+------------------+--------------+
| cudnn   | 5.1 | module load      | /cm/shared/apps/ | CUDA Deep    |
|         |     | cudnn/5.1        | cudnn/5.1        | Neural       |
|         |     |                  |                  | Network      |
+---------+-----+------------------+------------------+--------------+
| deal.ii | 8.4 | module load      | /rigel/opt/deal. | Library -    |
|         | .0  | deal.II/8.4.0    | II               | Scientific   |
+---------+-----+------------------+------------------+--------------+
| fftw2   | 2.1 | module load      | /cm/shared/apps/ | Library -    |
|         | .5  | fftw2/openmpi/gc | fftw/openmpi/gcc | Scientific   |
|         |     | c/64/double/2.1. | /64/2.1.5/double |              |
|         |     | 5                |                  |              |
+---------+-----+------------------+------------------+--------------+
| fftw3   | 3.3 | module load      | /cm/shared/apps/ | Library -    |
|         | .4  | fftw3/openmpi/gc | fftw/openmpi/gcc | Scientific   |
|         |     | c/64/3.3.4       | /64/3.3.4        |              |
+---------+-----+------------------+------------------+--------------+
| gcc     | 6.1 | module load      | /bin/gcc         | Compiler -   |
|         | .0  | gcc/6.1.0        |                  | C, C++       |
+---------+-----+------------------+------------------+--------------+
| gcc     | 4.8 | module load      | /bin/gcc         | Compiler -   |
|         | .5  | gcc/4.8.5        |                  | C, C++       |
+---------+-----+------------------+------------------+--------------+
| git     | 1.8 | /usr/bin/git     | /usr/bin/git     | Revision     |
|         | .3. |                  |                  | control      |
|         | 1   |                  |                  |              |
+---------+-----+------------------+------------------+--------------+
| gromacs | 5.1 | module load      | /rigel/opt/groma | Molecular    |
|         | .4  | gromacs/5.1.4-no | cs-5.1.4-nongpu  | Dynamics     |
|         |     | ngpu             |                  |              |
+---------+-----+------------------+------------------+--------------+
| gromacs | 5.1 | /rigel/opt/groma | /rigel/opt/groma | Molecular    |
| gpu     | .4  | cs-5.1.4         | cs-5.1.4         | Dynamics     |
+---------+-----+------------------+------------------+--------------+
| gsl     | 1.1 | /usr/bin/gsl-his | /usr/bin/gsl-his | Library -    |
|         | 5   | togram           | togram           | Scientific   |
+---------+-----+------------------+------------------+--------------+
| idl     | 8.4 | module load      | /rigel/opt/exeli | Library -    |
|         |     | idl/8.4          | s/idl84          | Scientific   |
+---------+-----+------------------+------------------+--------------+
| intel-p | 201 | module load      | /rigel/opt/paral | Intel        |
| arallel | 7   | intel-parallel-s | lel_studio_xe_20 | Compiler     |
| -studio |     | tudio/2017       | 17               |              |
+---------+-----+------------------+------------------+--------------+
| java    | 1.8 | module load      | /usr/bin/java    | Java         |
|         | .0  | java/1.8.0       |                  |              |
+---------+-----+------------------+------------------+--------------+
| matlab  | R20 | module load      | /rigel/opt/matla | Numerical    |
|         | 16b | matlab/2016b     | b/R2016b         | Computing    |
+---------+-----+------------------+------------------+--------------+
| netcdf- | 4.4 | module load      | /rigel/opt/netcd | Library -    |
| fortran | .4  | netcdf-fortran/4 | f-fortran-4.4.4  | Scientific   |
|         |     | .4.4             |                  |              |
+---------+-----+------------------+------------------+--------------+
| paralle | 1.8 | module load      | /rigel/opt/hdf5p | Library -    |
| l       | .18 | hdf5p/1.8.18     | -1.8.18          | Scientific   |
| hdf5    |     |                  |                  |              |
+---------+-----+------------------+------------------+--------------+
| R       | 3.3 | module load      | /rigel/opt/R-3.3 | Programming  |
|         | .2  | R/3.3.2          | .2               | language     |
+---------+-----+------------------+------------------+--------------+
| schrodi | 201 | module load      | /rigel/opt/schro | Chemical     |
| nger    | 6-3 | schrodinger/2016 | dinger/2016-3    | simulation   |
|         |     | -3               |                  |              |
+---------+-----+------------------+------------------+--------------+
| schrodi | 201 | module load      | /rigel/opt/schro | Chemical     |
| nger    | 6-4 | schrodinger/2016 | dinger/2016-4    | simulation   |
|         |     | -4               |                  |              |
+---------+-----+------------------+------------------+--------------+
| Stata   | 201 | module load      |                  | Stata        |
|         | 5   | stata            |                  |              |
+---------+-----+------------------+------------------+--------------+

`Habanero "Submitting Jobs" page <https://confluence.columbia.edu/confluence/display/rcs/Habanero+-+Submitting+Jobs/>`_
tells you how to submit batch jobs on Habanero.

`Habanero "Storage" page <https://confluence.columbia.edu/confluence/display/rcs/Habanero+-+Storage/>`_ tells you how large storage you can have on Habanero.
