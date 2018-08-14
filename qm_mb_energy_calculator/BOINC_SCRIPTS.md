# Instructions for all scripts relevent to BOINC

## General:
The workflow is broken down into 4 stages, each of which has a unique script

Many of the scripts require a settings file. A settings file contains additional information that is required by the script but not passed in via arguments. By convention, settings files should end in ".ini". The required fields for a settings file to get the below scripts to work is as such:

```
[files]
log_path = \some\path\where\logs\should\go

[molcule]
fragments = 1,2,3 # comma delimerated list of how many atoms are in each fragment
charges = 0,-1,0 # comma delimerated list of the charges of each fragment
spins = 1,2,1 # comma delimerated list of the spin multiplicity of each fragment
names = H,CN-,H2O # comma delimerated list of the names of each fragment
tag = someTag # you can specify energies by their tag when generating a training set

[energy_calculator]
method = HF # method to use for calculations
basis = STO-3G # basis to use for calculations
cp = True # whether to use counterpoise correction for calculations

[psi4]
memory = 500MB # amount of memory to use for calculations. can be specified in MB, GB, etc
num_threads = 6 # maximum number of threads, will be reduced to the maximum number of threads a machine can support
```

## Creating a database:
A database is created with database_initializer.py:

database_initializer.py creates a database from a set of configurations, or updates an existing database to contain new configurations

```
python database_initializer.py <settings file> <database name> <config files>
```

settings file: apropriate ".ini" file as discussed in the "General" section above

database name: filepath to where the database should be created

config files: either a directory containing many ".xyz" files, or the path to a signle ".xyz" file. All ".xyz" files will be initialized into the database, other files will be ignored. Any file ending in ".opt.xyz" will be considered an optimized geometry and marked in the database as such.

## Creating a job
Jobs are created with database_job_maker.py:

database_job_maker.py makes one job for every uncalculated energy in a database

```
python database_job_maker.py <settings file> <database name> <job directory>
```

settings file: apropriate ".ini" file as discussed in the "General" section above

database name: the database file

job directory: the jobs will be created in this directory

## Running the jobs
After running database_job_maker.py, the supplied job directory will be filled with .py Job files. Each of these files is a stand-alone python script that can be run with the "python" command like any other python script.

Each Job will generate a .out file and a .log file. These are used to enter information back into the database.

## Reading a job
Jobs are read by database_job_reader.py:

database_job_reader.py takes in a job's .log and .out file and enters the information into a database

```
python database_job_reader.py <database name> <job out path> <job log path>
```

database name: the database file

job out path: the job's .out file

job log path: the job's .log file
