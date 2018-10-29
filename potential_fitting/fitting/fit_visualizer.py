
from potential_fitting.database import Database

import subprocess 

import matplotlib.pyplot as plt 

import numpy as np

def compare_energies(file_path_TTM, file_path_TTM_params, file_path_MB, file_path_MB_params, db_name, molecule_name, method, basis, cp, tag):
	ttm = []
	mb = []

	with Database() as database:

		energy_molecule_pairs = database.get_energies(molecule_name, method, basis, cp, tag)

		molecules = [i[0] for i in energy_molecule_pairs]

		calc = [j[1] for j in energy_molecule_pairs]

		for m in molecules:
			molecule_xyz = m.to_xyz()

			with open("file1.txt", "w") as f:
				f.write(m.get_num_atoms())
				f.write('\n')
				f.write("######")
				f.write(molecule_xyz)

			result_ttm = subprocess.run([file_path_TTM, file_path_TTM_params, "file1.txt"], stdout=subprocess.PIPE)

			ttm += [float(result_ttm.stdout.split[2])]

			result_mb = subprocess.run([file_path_MB, file_path_MB_params, "file1.txt"], stdout=subprocess.PIPE)

			mb += [float(result_mb.stdout.split[2])]



		ttm_array = np.asarray(ttm)
		calc_array = np.asarray(calc)
		mb_array = np.asarray(mb)


		plt.scatter(calc_array, ttm_array)
		plt.scatter(calc_array, mb_array)
		plt.show()

		plt.clf()

		plt.scatter(calc_array, ttm_array - calc_array)
		plt.scatter(calc_array, mb_array - calc_array)
		plt.show()
















