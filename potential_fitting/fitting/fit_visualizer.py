
from potential_fitting.database import Database

import subprocess 

import matplotlib.pyplot as plt 

from potential_fitting.utils import constants

import numpy as np

def compare_energies(file_path_TTM, file_path_TTM_params, file_path_MB, file_path_MB_params, db_name, molecule_name, method, basis, cp, tag):
    ttm = []
    mb = []

    with Database(db_name) as database:

        energy_molecule_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))
    
        molecules = [i[0] for i in energy_molecule_pairs]

        calc = [j[1][2] for j in energy_molecule_pairs]
        
        min_calc = min(calc)

        calc = [(i - min_calc) * constants.au_to_kcal for i in calc]

        for m in molecules:
            molecule_xyz = m.to_xyz()

            with open("file1.txt", "w") as f:
                f.write(str(m.get_num_atoms()))
                f.write('\n')
                f.write("######\n")
                f.write(molecule_xyz)

            result_ttm = subprocess.run([file_path_TTM, file_path_TTM_params, "file1.txt"], stdout=subprocess.PIPE)

            ttm += [float(result_ttm.stdout.split()[2])]

            result_mb = subprocess.run([file_path_MB, file_path_MB_params, "file1.txt"], stdout=subprocess.PIPE)

            mb += [float(result_mb.stdout.split()[2])]



        ttm_array = np.asarray(ttm)
        calc_array = np.asarray(calc)
        mb_array = np.asarray(mb)

        plt.figure(1)
        ttm_vs_calc = plt.scatter(calc_array, ttm_array)
        mb_vs_calc = plt.scatter(calc_array, mb_array)
        plt.plot(calc_array, calc_array, 'r--')
        plt.legend((ttm_vs_calc, mb_vs_calc), ('TTM', 'MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Fitted Energy [Kcal/mol]")
        #plt.show()

        #plt.clf()
        plt.figure(2)
        error_ttm = plt.scatter(calc_array, ttm_array - calc_array)
        error_mb = plt.scatter(calc_array, mb_array - calc_array)
        plt.plot(calc_array, [0 for i in range(len(calc_array))], 'r--')
        plt.legend((error_ttm, error_mb), ('Error TTM', 'Error MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Error [Kcal/mol]")

        plt.show()
















