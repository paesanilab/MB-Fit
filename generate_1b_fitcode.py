import sys
import os
import argparse
import subprocess

def execute_with_output(cmd, output):
    f = open(output, "w")
    subprocess.check_call(cmd, stdout=f)
    f.close()

def execute(cmd):
    subprocess.check_call(cmd)
    #popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    #for stdout_line in iter(popen.stdout.readline, ""):
    #    yield stdout_line 
    #popen.stdout.close()
    #return_code = popen.wait()
    #if return_code:
    #    raise subprocess.CalledProcessError(return_code, cmd)

parser = argparse.ArgumentParser(description='Run mbpol polynomial generation and fitting')

parser.add_argument('configuration_file_path',
                    help='path to the configuration file')
args = parser.parse_args()



import configparser

config = configparser.ConfigParser()
config.read(args.configuration_file_path)


number_of_bodies = len(config.get("common", "molecule").split("_"))


# # Polynomial generation

order = config.getint("common", "polynomial_order")


input_file = config.get("common", "molecule")

execute(["generate_input_poly.py", input_file])

os.makedirs("polynomial_generation", exist_ok=True)
os.chdir("polynomial_generation")

poly_output_file = "poly.log"
execute_with_output(["poly-gen_mb-nrg.pl", str(order), '../' + input_file], poly_output_file)

for f in ["poly-nogrd", "poly-grd"]:
    maple_file = f + ".maple"
    execute(['maple', maple_file])
    c_file = f + ".c"
    cpp_file = f + ".cpp"
    with open(c_file, 'r') as c_file_handle, open(cpp_file, 'w') as cpp_file_handle:
        subprocess.check_call(['clean-maple-c.pl'], stdin=c_file_handle, stdout=cpp_file_handle)

# # Fitting

os.chdir("..")
os.makedirs("fitting", exist_ok=True)
os.chdir("fitting")


# ### Setup fitting and compile

execute(['prepare_1b_fitting_code.sh', '../'+input_file, '../polynomial_generation/', '../config.ini'])

# Commented, unnecessary
#execute(["./fit-1b", config.get("fitting", "training_set_file")])
