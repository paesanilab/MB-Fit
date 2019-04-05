import psi4
import subprocess, os

"""
This is not a runnable python file.

It is a template that database_job_maker.py uses to make psi4 jobs
"""

def execute_job():
    whole_molecule = "{whole_molecule}"
    molecule = "{molecule}"
    frag_indices = "{frag_indices}"
    method = "{method}"
    basis = "{basis}"
    cp = "{cp}"
    number_of_threads = {num_threads}
    memory = "{memory}"

    try:
        max_threads = int(subprocess.check_output(["grep", "-c", "cores", "/proc/cpuinfo"]))
        print("Maximum threads: {format}".format(max_threads))
    except:
        print("Error detecting number of cores. \n Maxium threads: 1")
        max_threads = 1

    if number_of_threads > max_threads:
        print("Input number of threads ({format}) greater than max threads ({format}), limiting number of threads to max threads".format(number_of_threads, max_threads))
        number_of_threads = max_threads

    print("Running Job")
    print("Molcule: {format}".format(molecule))
    print("Method: {format}".format(method))
    print("Basis: {format}".format(basis))
    print("Threads {format}".format(number_of_threads))
    print("Memory: {format}".format(memory))

    i = 1

    job_dir =  "job_{format}".format(i)

    while os.path.isfile(job_dir):

        i += 1

        job_dir = "job_{format}".format(i)

    os.mkdir(job_dir)
    output_file = job_dir + "/output.dat"
    log_file = job_dir + "/output.log"

    psi4.core.set_output_file(log_file, False)
    psi4.set_memory(memory)
    psi4.geometry(molecule)
    psi4.set_num_threads(number_of_threads)

    try:
        energy = psi4.energy("{format}/{format}".format(method, basis))
        print("Energy: {format}".format(energy))
        success = True
    except ValueError:
        success = False
        print("Iterations failed to Converge")


    with open(output_file, "w") as out_file:
        out_file.writelines(["Molecule: {format}".format(whole_molecule),
                             "\nMethod: {format}".format(method),
                             "\nBasis: {format}".format(basis),
                             "\nCp: {format}".format(cp),
                             "\nfrag_indices: {format}".format(frag_indices)
                             ])

        if success:
            out_file.write("\nSuccess: {format}".format(energy))
        else:
            out_file.write("\nFailure")

if __name__ == "__main__":
    execute_job()
