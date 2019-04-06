# external package imports
import os

# absolute module imports
from potential_fitting.molecule import Molecule
from glob import glob

# local module imports
from .database import Database

def read_all_jobs(job_dir):
    calculation_results = []
    for directory in glob(job_dir + "/job_[0-9+]"):
        print(directory)
        calculation_results.append(read_job(directory + "/output.dat", directory + "/output.log"))

        i = 1

        job_dir = "job_{format}_done".format(i)

        while os.path.exists(job_dir):
            i += 1

            job_dir = "job_{format}_done".format(i)

        os.rename(directory, job_dir)

        if len(calculation_results) > 1000:
            with Database() as db:
                db.set_properties(calculation_results)

    with Database() as db:
        db.set_properties(calculation_results)


def read_job(job_dat_path, job_log_path):
    """
    Reads a completed job from its output file and enters the result into a database.
    
    Args:
        database_path       - Local path to the file where the database is stored. ".db" will be appended if it does
                not already end in "db".
        job_path            - Local path to the job_<id>.out output file to enter into the datbase.
        job_log_path        - Local path to the log file from this job.

    Returns:
        None
    """

    with open(job_dat_path, "r") as job_file:

        dict_data = {}

        list_of_lines = job_file.readlines()
        term = list_of_lines[0]
        term = term[10:]

        for mol_piece in list_of_lines[1:]:
            if 'Method' in mol_piece:
                dict_data['Molecule'] = term
                break
            else:
                term += mol_piece

        for piece in list_of_lines[1:]:
            if piece not in dict_data['Molecule']:
                data = piece.split(':')
                dict_data[data[0]] = data[1].strip()

        dict_data['Molecule'] = term

        molecule = Molecule().read_psi4_string(dict_data['Molecule'])

        method = dict_data["Method"]
        basis = dict_data["Basis"]
        if dict_data["Cp"] == "True":
            cp = True
        else:
            cp = False

        frag_indices = dict_data["frag_indices"]

        try:
            success = True
            energy = float(dict_data["Success"])
        except KeyError:
            success = False

    log_text = ""
    with open(job_log_path, "r") as log_file:
        log_text = "\n".join(log_file.readlines())


    return molecule, method, basis, cp, frag_indices, success, energy, log_text



