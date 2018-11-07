
from potential_fitting.database import Database

import subprocess 

import matplotlib.pyplot as plt 

from potential_fitting.utils import constants

import numpy as np

def rmsd(error_array):
    rmsd_calculated = np.sqrt(np.mean(error_array**2))
    return rmsd_calculated
    

def compare_energies(file_path_TTM, file_path_TTM_params, file_path_MB, file_path_MB_params, db_name, monomer1_name, monomer2_name, method, basis, cp, tag, threshold = 50):
    ttm = []
    mb = []

    molecule_name = "-".join(sorted([monomer1_name, monomer2_name]))

    with Database(db_name) as database:

        energy_molecule_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))

        monomer1_opt = list(database.get_energies(monomer1_name, method, basis, cp, tag))[0][1][0]
        monomer2_opt = list(database.get_energies(monomer2_name, method, basis, cp, tag))[0][1][0]
    
        molecules = [i[0] for i in energy_molecule_pairs]

        calc = [j[1][2] - j[1][1] - j[1][0] for j in energy_molecule_pairs]

        binding_energies = [interaction - 
                ((energies[1][0] if len(energies[1]) == 3 else energies[1][3]) - monomer1_opt) - 
                ((energies[1][1] if len(energies[1]) == 3 else energies[1][4]) - monomer2_opt) 
                for interaction, energies in zip(calc, energy_molecule_pairs)]

        binding_energies = [i * constants.au_to_kcal for i in binding_energies]

        calc = [i * constants.au_to_kcal for i in calc]

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
            if binding_energies[i] >= threshold:
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

        ttm_vs_calc_above = plt.scatter(calc_above_array, ttm_above_array, c = '#00074c', s = 5, alpha = 0.5)
        ttm_vs_calc_below = plt.scatter(calc_below_array, ttm_below_array, c = '#0017ff', s = 5, alpha = 0.5)

        

        mb_vs_calc_above = plt.scatter(calc_above_array, mb_above_array, c = '#db2000', s = 5, alpha = 0.5 )
        mb_vs_calc_below = plt.scatter(calc_below_array, mb_below_array, c = '#5b0d00', s = 5, alpha = 0.5)

        plt.plot(calc_array, calc_array, c = 'orange', alpha = 0.5)
        plt.legend((ttm_vs_calc_above, mb_vs_calc_above), ('TTM', 'MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Fitted Energy [Kcal/mol]")
        

        plt.figure(2)

        error_ttm_above = plt.scatter(calc_above_array, ttm_above_array - calc_above_array, c = '#00074c', s = 5, alpha = 0.5)
        error_ttm_below = plt.scatter(calc_below_array, ttm_below_array - calc_below_array, c = '#0017ff', s = 5, alpha = 0.5)

        error_mb_above = plt.scatter(calc_above_array, mb_above_array - calc_above_array, c = '#db2000', s = 5, alpha = 0.5)
        error_mb_below = plt.scatter(calc_below_array, mb_below_array - calc_below_array, c = '#5b0d00', s = 5, alpha = 0.5)

        rmsd_ttm = rmsd(ttm_array - calc_array)
        rmsd_mb = rmsd(mb_array - calc_array)

        rmsd_positive_ttm = np.array([rmsd_ttm for i in calc_array])
        rmsd_negative_ttm = np.array([-rmsd_ttm for i in calc_array])

        rmsd_positive_mb = np.array([rmsd_mb for i in calc_array])
        rmsd_negative_mb = np.array([-rmsd_mb for i in calc_array])

        plt.plot(calc_array, [0 for i in calc_array], c = 'orange', alpha = 0.5)
        plt.plot(calc_array, rmsd_positive_ttm, c= 'blue')
        plt.plot(calc_array, rmsd_negative_ttm, c = 'blue')

        plt.plot(calc_array, rmsd_positive_mb, c = 'red')
        plt.plot(calc_array, rmsd_negative_mb, c = 'red')

        plt.legend((error_ttm_above, error_mb_above), ('Error TTM', 'Error MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Error [Kcal/mol]")


        plt.figure(3)
        
        ttm_vs_calc_below = plt.scatter(calc_below_array, ttm_below_array, c = '#0017ff', s = 5, alpha = 0.5)

        mb_vs_calc_below = plt.scatter(calc_below_array, mb_below_array, c = '#5b0d00', s = 5, alpha = 0.5)
        
        plt.plot(calc_below_array, calc_below_array, c = 'orange', alpha = 0.5)


        plt.legend((ttm_vs_calc_below, mb_vs_calc_below), ('TTM', 'MB'))

        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Fitted Energy [Kcal/mol]")

        plt.figure(4)

        
        error_ttm_below = plt.scatter(calc_below_array, ttm_below_array - calc_below_array, c = '#0017ff', s = 5, alpha = 0.5)

        error_mb_below = plt.scatter(calc_below_array, mb_below_array - calc_below_array, c = '#5b0d00', s = 5, alpha = 0.5)



        rmsd_ttm = rmsd(ttm_below_array - calc_below_array)
        rmsd_mb = rmsd(mb_below_array - calc_below_array)

        rmsd_positive_ttm = np.array([rmsd_ttm for i in calc_below_array])
        rmsd_negative_ttm = np.array([-rmsd_ttm for i in calc_below_array])

        rmsd_positive_mb = np.array([rmsd_mb for i in calc_below_array])
        rmsd_negative_mb = np.array([-rmsd_mb for i in calc_below_array])

        plt.plot(calc_below_array, rmsd_positive_ttm, c= 'blue')
        plt.plot(calc_below_array, rmsd_negative_ttm, c = 'blue')

        plt.plot(calc_below_array, rmsd_positive_mb, c = 'red')
        plt.plot(calc_below_array, rmsd_negative_mb, c = 'red')



        plt.plot(calc_below_array, [0 for i in calc_below_array], c = 'orange', alpha = 0.5)

        
        plt.legend((error_ttm_above, error_mb_above), ('Error TTM', 'Error MB'))
        plt.xlabel("Ref. Energy [Kcal/mol]")
        plt.ylabel("Error [Kcal/mol]")

        plt.show()
















