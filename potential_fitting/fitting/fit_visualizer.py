
from potential_fitting.database import Database

import subprocess 

import matplotlib.pyplot as plt 

from potential_fitting.utils import constants

import numpy as np

class Dataset():

    colors = [
        [[[1, 0, 0]], [[1, 0.5, 0.5]]],
        [[[0, 1, 0]], [[0.5, 1, 0.5]]],
        [[[0, 0, 1]], [[0.5, 0.5, 1]]]
    ]

    def __init__(self, calc_energies, fit_energies, method):
        self.calc_energies = calc_energies
        self.fit_energies = fit_energies
        self.method = method

    def split_at_threshold(self, threshold):
        raise NotImplementedError

class Dataset_1b(Dataset):

    def split_at_threshold(self, threshold):
        low_calc = []
        high_calc = []

        low_fit = []
        high_fit = []

        for calc_energy, fit_energy in zip(self.calc_energies, self.fit_energies):
            if calc_energy < threshold:
                low_calc.append(calc_energy)
                low_fit.append(fit_energy)
            else:
                high_calc.append(calc_energy)
                high_fit.append(fit_energy)

        return Dataset_1b(low_calc, low_fit, self.method + " below {} kcal/mol".format(threshold)), Dataset_1b(high_calc, high_fit, self.method + " above {} kcal/mol".format(threshold))


class Dataset_2b(Dataset):

    def __init__(self, calc_energies, fit_energies, method, binding_energies):
        super(Dataset_2b, self).__init__(calc_energies, fit_energies, method)

        self.binding_energies = binding_energies

    def split_at_threshold(self, threshold):
        low_calc = []
        high_calc = []

        low_fit = []
        high_fit = []

        low_bind = []
        high_bind = []

        for calc_energy, fit_energy, binding_energy in zip(self.calc_energies, self.fit_energies, self.binding_energies):
            if calc_energy < threshold:
                low_calc.append(calc_energy)
                low_fit.append(fit_energy)
                low_bind.append(binding_energy)
            else:
                high_calc.append(calc_energy)
                high_fit.append(fit_energy)
                high_bind.append(binding_energy)

        return Dataset(low_calc, low_fit, self.method + " below {} kcal/mol".format(threshold), low_bind), Dataset(high_calc, high_fit, self.method + " above {} kcal/mol".format(threshold), high_bind)

def make_1b_graphs(file_path_MB, file_path_MB_params, db_name, molecule_name, method, basis, cp, tag,
        low_threshold = 50):
    
    mb = []
    with Database(db_name) as database:

        energy_molecule_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))

        # getting the monomer optimized energies 
        opt_energy = list(database.get_energies(molecule_name, method, basis, cp, tag))[0][1][0]

    # getting the molecules
    molecules = [i[0] for i in energy_molecule_pairs]

    # calculating the required energy from the energy-molecule pairs
    calc = [i[1][0] - opt_energy for i in energy_molecule_pairs]

    calc = [i * constants.au_to_kcal for i in calc]

    for m in molecules:
        #writing to xyz file
        molecule_xyz = m.to_xyz()

        with open("file1.txt", "w") as f:

            f.write(str(m.get_num_atoms()))
            f.write('\n')
            f.write("######\n")
            f.write(molecule_xyz)

        #getting the mb data
        result_mb = subprocess.run([file_path_MB, file_path_MB_params, "file1.txt"], stdout=subprocess.PIPE)
        #adding mb_data to a list
        mb += [float(result_mb.stdout.split()[2])]

    make_graphs(Dataset_1b(calc, mb, "{}/{}/{}".format(method, basis, cp)), low_threshold = low_threshold)

def make_2b_graphs(file_path_ttm, file_path_ttm_params, file_path_MB, file_path_MB_params, db_name, monomer1_name,
        monomer2_name, method, basis, cp, tag, low_threshold = 50):

    ttm = []
    mb = []

    #creating a dimer
    molecule_name = "-".join(sorted([monomer1_name, monomer2_name]))

    with Database(db_name) as database:

        energy_molecule_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))

        #getting the monomer optimized energies 
        monomer1_opt = list(database.get_energies(monomer1_name, method, basis, cp, tag))[0][1][0]
        monomer2_opt = list(database.get_energies(monomer2_name, method, basis, cp, tag))[0][1][0]

        #getting the molecules
        molecules = [i[0] for i in energy_molecule_pairs]

         #calculating the required interaction energy from the energy-molecule pairs
        calc = [j[1][2] - j[1][1] - j[1][0] for j in energy_molecule_pairs]


        #getting the binding energies
        binding_energies = [interaction - 
                ((energies[1][0] if len(energies[1]) == 3 else energies[1][3]) - monomer1_opt) - 
                ((energies[1][1] if len(energies[1]) == 3 else energies[1][4]) - monomer2_opt) 
                for interaction, energies in zip(calc, energy_molecule_pairs)]

        binding_energies = [i * constants.au_to_kcal for i in binding_energies]

        calc = [i * constants.au_to_kcal for i in calc]

        for m in molecules:
            #writing to xyz file
            molecule_xyz = m.to_xyz()

            with open("file1.txt", "w") as f:

                f.write(str(m.get_num_atoms()))
                f.write('\n')
                f.write("######\n")
                f.write(molecule_xyz)

            result_ttm = subprocess.run([file_path_TTM, file_path_TTM_params, "file1.txt"], stdout=subprocess.PIPE)

            #getting the ttm_data and adding to a list
            ttm += [float(result_ttm.stdout.split()[2])]

            #getting the mb data
            result_mb = subprocess.run([file_path_MB, file_path_MB_params, "file1.txt"], stdout=subprocess.PIPE)
            #adding mb_data to a list
            mb += [float(result_mb.stdout.split()[2])]

    make_graphs(Dataset_2b(calc, ttm, "ttm", binding), Dataset_2b(calc, mb, "{}/{}/{}".format(method, basis, cp), binding), low_threshold = low_threshold)


def make_graphs(*datasets, low_threshold = 50):
    make_energy_graph(1, *datasets, low_threshold = low_threshold)
    make_error_graph(2, *datasets, low_threshold = low_threshold)

def make_energy_graph(figure_num, *datasets, low_threshold = 50):
    # make the graph featuring all information divided into low and high energy datasets

    # Set figure number
    plt.figure(figure_num)

    above_plots = []
    below_plots = []

    # plot each dataset
    for index, dataset in enumerate(datasets):

        low_dataset, high_dataset = dataset.split_at_threshold(low_threshold)

        below_plots.append(plt.scatter(low_dataset.calc_energies, low_dataset.fit_energies, c = Dataset.colors[index][0], s = 5, alpha = 0.5))
        above_plots.append(plt.scatter(high_dataset.calc_energies, high_dataset.fit_energies, c = Dataset.colors[index][1], s = 5, alpha = 0.5))

    # plotting an idealized prediction using color codes for TTM fit
    # NOT IDEAL, should just plot y=x constrained to the graph
    plt.plot(datasets[0].calc_energies, datasets[0].calc_energies, c = 'orange', alpha = 0.5)

    #Adding a legend
    plt.legend([plot for plot in below_plots], [dataset.method for dataset in datasets])

    #Adding axes titles
    plt.xlabel("Ref. Energy [Kcal/mol]")
    plt.ylabel("Fitted Energy [Kcal/mol]")

    # make the graph of the error featuring all info divided into low and high energy datasets

def make_error_graph(figure_num, *datasets, low_threshold = 50):
    # make the graph featuring all error information high and low

    plt.figure(figure_num)

    above_plots = []
    below_plots = []

    for index, dataset in enumerate(datasets):

        low_dataset, high_dataset = dataset.split_at_threshold(low_threshold)

        below_plots.append(plt.scatter(low_dataset.calc_energies,
                [fit - calc for fit, calc in zip(low_dataset.fit_energies, low_dataset.calc_energies)],
                c = Dataset.colors[index][0], s = 5, alpha = 0.5))
        below_plots.append(plt.scatter(high_dataset.calc_energies,
                [fit - calc for fit, calc in zip(high_dataset.fit_energies, high_dataset.calc_energies)],
                c = Dataset.colors[index][1], s = 5, alpha = 0.5))
    
    # plotting an idealized prediction using color codes for TTM fit
    # NOT IDEAL, should just plot y=x constrained to the graph
    plt.plot(datasets[0].calc_energies, [0 for calc_energy in datasets[0].calc_energies], c = 'orange', alpha = 0.5)

    #Adding a legend
    plt.legend([plot for plot in below_plots], [dataset.method for dataset in datasets])

    #Adding axes titles
    plt.xlabel("Ref. Energy [Kcal/mol]")
    plt.ylabel("Fitted Energy - Ref. Energy [Kcal/mol]")

    # make the graph of the error featuring all info divided into low and high energy datasets


def rmsd(error_array):

    '''
    Calculates and returns the root-mean-square-distance, given a numpy array

    Args:
        error_array: The array whose rmsd needs to be calculated.
    '''

    rmsd_calculated = np.sqrt(np.mean(error_array**2))

    return rmsd_calculated
