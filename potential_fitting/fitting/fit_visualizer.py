
from potential_fitting.database import Database

import subprocess 

import matplotlib.pyplot as plt 

from potential_fitting.utils import constants

import numpy as np

def compare_energies(file_path_TTM, file_path_TTM_params, file_path_MB, file_path_MB_params, db_name, molecule_name, method, basis, cp, tag, threshold = 50):
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


        calc_below = []
        calc_above = []

        ttm_below = []
        ttm_above = []

        mb_below = []
        mb_above = []


        for i in range(len(calc)):
            if calc[i] >= 50:
                calc_above += [calc[i]]
                ttm_above += [ttm[i]]
                mb_above += [mb[i]]
            else:
                calc_below += [calc[i]]
                ttm_below += [ttm[i]]
                mb_below += [mb[i]]

        ttm_array = np.asarray(ttm)
        calc_array = np.asarray(calc)
        mb_array = np.asarray(mb)

        calc_above_array = np.asarray(calc_above)
        calc_below_array = np.asarray(calc_below)

        ttm_above_array = np.asarray(ttm_above)
        ttm_below_array = np.asarray(ttm_below)

        mb_above_array = np.asarray(mb_above)
        mb_below_array = np.asarray(mb_below)



        plt.figure(1)

        ttm_vs_calc_above = plt.scatter(calc_above_array, ttm_above_array, c = '#00008B' )
        ttm_vs_calc_below = plt.scatter(calc_below_array, ttm_below_array, c = '#00BFFF')

        mb_vs_calc_above = plt.scatter(calc_above_array, mb_above_array, c = '#008000' )
        mb_vs_calc_below = plt.scatter(calc_below_array, mb_below_array, c = '##00FA9A')

        plt.plot(calc_array, calc_array, 'r--')
        plt.legend((ttm_vs_calc_above, mb_vs_calc_above), ('TTM', 'MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Fitted Energy [Kcal/mol]")
        

        plt.figure(2)

        error_ttm_above = plt.scatter(calc_above_array, ttm_above_array - calc_above_array, c = '#00008B')
        error_ttm_below = plt.scatter(calc_below_array, ttm_below_array - calc_below_array, c = '#00BFFF')

        error_mb_above = plt.scatter(calc_above_array, mb_above_array - calc_above_array, c = '#008000')
        error_mb_below = plt.scatter(calc_below_array, mb_below_array - calc_below_array, c = '##00FA9A')


        plt.plot(calc_array, [0 for i in range(len(calc_array))], 'r--')
        plt.legend((error_ttm_above, error_mb_above), ('Error TTM', 'Error MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Error [Kcal/mol]")

        plt.show()
















