import psi4
import subprocess

"""
This is not a runnable python file.

It is a template that database_job_maker.py uses to make psi4 jobs
"""

def execute_job():
    job_id = {job_id}
    molecule = "{molecule}"
    method = "{method}"
    basis = "{basis}"
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

    psi4.core.set_output_file("job_{format}.log".format(job_id), False)
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

    with open("job_{format}.out".format(job_id), "w") as out_file:
        out_file.write("Job: {format}\n".format(job_id))
        if success:
            out_file.write("Energy: {format}".format(energy))
        else:
            out_file.write("Failure")

if __name__ == "__main__":
    execute_job()
