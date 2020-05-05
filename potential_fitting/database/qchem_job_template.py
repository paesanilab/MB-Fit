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
    job_hash = "{job_hash}"
    arguments = {arguments}

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

    i = 8
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "job_{format}".format(job_hash[:i]))

    while os.path.exists(job_dir):
        i += 1

        job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "job_{format}".format(job_hash[:i]))

    os.mkdir(job_dir)
    input_path = job_dir + "/qchem_input.log"
    qchem_stdout_stderr_log_path = job_dir + "/qchem_log.log"
    output_path = job_dir + "/output.ini"
    log_path = job_dir + "/output.log"

    # qchem specific stuff

    with open(input_path, "w") as input_file:
        input_file.write("$molecule\n")
        input_file.write("{{}} {{}}\n".format(total_charge, total_spin))
        input_file.write("{{}}\n".format(molecule))
        input_file.write("$end\n")

        input_file.write("$rem\n")
        input_file.write("jobtype sp\n")
        input_file.write("method {{}}\n".format(method))
        input_file.write("basis {{}}\n".format(basis))

        for key, value in arguments.items():
            input_file.write(str(key) + " " + str(value))

        input_file.write("$end\n")


    subprocess.run(["touch", qchem_stdout_stderr_log_path])
    with open(qchem_stdout_stderr_log_path, "r") as qchem_stdout_stderr_log_file:
        subprocess.run(["qchem", "-nt", str(number_of_threads), input_path, log_path], stdout=qchem_stdout_stderr_log_file, stderr=qchem_stdout_stderr_log_file)

    energy = None
    success = False
    with open(log_path) as out_file:
        for line in out_file:
            if line.find("Total energy in the final basis set = ") != -1:
                energy = float(line[line.find("Total energy in the final basis set = ") + 39:])
                success = True
                print("Energy: {format}".format(energy))
                break

    # end qchem specific stuff

    with open(output_path, "w") as out_file:
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
        out_file.write("job_hash = {format}\n".format(job_hash))

        if success:
            out_file.write("energy = {format}\n".format(energy))


if __name__ == "__main__":
    execute_job()
