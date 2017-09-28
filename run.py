import sys
import os
import argparse
import subprocess

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

execute(["poly-gen_mb-nrg.pl", str(order), '../' + input_file])

for f in ["poly-nogrd", "poly-grd"]:
    maple_file = f + ".maple"
    execute(['maple', maple_file])
    c_file = f + ".c"
    cpp_file = f + ".cpp"
    execute(['clean-maple-c.pl', c_file, cpp_file])

# # Fitting

os.chdir("..")
os.makedirs("fitting", exist_ok=True)
os.chdir("fitting")


# ### Setup fitting and compile

execute(['prepare_1b_fitting_code.sh', '../'+input_file, '../polynomial_generation/', "/home/azonca/Paesani/potential_fitting/fitting/1B/get_codes/template", '../config.ini'])

execute(["./fit-1b", config.get("fitting", "training_set_file")])
