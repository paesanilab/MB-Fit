import sys
import configparser

def generate_poly(settings):
    
    # create configparser
    config = configparser.SafeConfigParser(allow_no_value=False)
    config.read(settings)

    # order of polynomial to generate
    poly_order = config["poly_generator"]["order"]
    input_file = config["files"]["poly_in_path"]

    # open output files for polynomials, will automatically be closed by with open as syntax
    with open(config["files"]["poly_path"] + "/vars.cpp", "w") as out_vars, open(config["files"]["poly_path"] + "/poly-direct.cpp", "w") as out_cpp, open(config["files"]["poly_path"] + "/poly-nogrd.maple", "w") as out_maple_nogrd, open(config["files"]["poly_path"] + "/poly-grd.maple", "w") as out_maple_grd, open(input_file, "r") as in_poly:
        
        # read the add_molecule definitions to get a list of the fragments
        fragments = parse_fragments(in_poly)

        # Parse the molecule names to determine and name the atoms in the system

        atom_names = []
        number_per_fragment = []
        index_this_fragment = []
        fragname = 'a'
        
        for fragment in fragments:
            pass

def parse_fragments(input_file):
    fragments = []
    line = input_file.readline()
    while line[:5] == "add_m":
        fragments.append(line[line.index("['") + 2:line.index("']")])
        line = input_file.readline()

    return fragments
   
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_poly.py settings.ini")
        sys.exit(1)
    generate_poly(sys.argv[1])
