import psi4
import subprocess, os

"""
This is not a runnable python file.

It is a template that database_job_maker.py uses to make psi4 jobs
"""

def execute_job():
    whole_molecule = "{whole_molecule}"
    charges = "{charges}"
    spins = "{spins}"
    symmetries = "{symmetries}"
    SMILES = "{SMILES}"
    names = "{names}"
    atom_counts = "{atom_counts}"
    total_atoms = "{total_atoms}"
    molecule = "{molecule}"
    frag_indices = "{frag_indices}"
    method = "{method}"
    basis = "{basis}"
    cp = "{cp}"
    use_cp = "{use_cp}"
    number_of_threads = {num_threads}
    memory = "{memory}"
    total_charge = "{total_charge}"
    total_spin = "{total_spin}"

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

    while os.path.exists(job_dir):

        i += 1

        job_dir = "job_{format}".format(i)

    os.mkdir(job_dir)
    output_file = job_dir + "/output.ini"
    log_file = job_dir + "/output.log"

    psi4.core.set_output_file(log_file, False)
    psi4.set_memory(memory)
    psi4_input_geo = "\n" + str(total_charge) + " " + str(total_spin) + "\n" + molecule
    psi4.geometry(psi4_input_geo)
    psi4.set_num_threads(number_of_threads)

    try:
        energy = psi4.energy("{format}/{format}".format(method, basis))
        print("Energy: {format}".format(energy))
        success = True
    except ValueError:
        success = False
        print("Iterations failed to Converge")


    with open(output_file, "w") as out_file:
        out_file.write("[molecule]\n")
        out_file.write("xyz = {format}\n\n{format}".format(total_atoms, whole_molecule).replace("\n", "\n ") + "\n")
        out_file.write("atom_counts = {format}\n".format(atom_counts))
        out_file.write("charges = {format}\n".format(charges))
        out_file.write("spins = {format}\n".format(spins))
        out_file.write("symmetries = {format}\n".format(symmetries).replace("'", ""))
        out_file.write("SMILES = {format}\n".format(SMILES).replace("'", ""))
        out_file.write("names = {format}\n".format(names).replace("'",""))
        out_file.write("method = {format}\n".format(method))
        out_file.write("basis = {format}\n".format(basis))
        out_file.write("cp = {format}\n".format(cp))
        out_file.write("use_cp = {format}\n".format(use_cp))
        out_file.write("frag_indices = {format}\n".format(frag_indices))

        if success:
            out_file.write("energy = {format}\n".format(energy))


if __name__ == "__main__":
    execute_job()
