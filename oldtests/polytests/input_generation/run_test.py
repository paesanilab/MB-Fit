import sys
import potential_fitting

potential_fitting.generate_poly_input("settings.ini", sys.argv[1].split(".")[0], sys.argv[1])
