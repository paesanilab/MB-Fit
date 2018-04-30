"""
A module for computing the energy of a molecule with or without
many-body decomposition.
"""

"""
LIBRARY IMPORTS
"""

import itertools
import calculator
from math import factorial

"""
Returns a 3d array of indicies to be passed into calc_energy.
Its hard to explain what this does, so I'll just give examples:    


if mbdecomp is true, returns an array of every possible combinations
in the list, sorted by size.
    Example: index_list = [0, 1, 2]
            return value = [
                [[0], [1], [2]],
                [[0, 1], [0, 2], [1, 2]],
                [[0, 1, 2]]
            ]

if mbdecomp is false:
    Example: index_list = [0, 1, 2]
            return value = [[[0, 1, 2]]]
"""
def build_frag_indices(index_list, mbdecomp):
    # used to hold array of combinations of index_list
    combinations_arr = []
    
    # if mbdecomp is False, just return a copy of index_list
    if mbdecomp == False:
        combinations_arr.append([index_list[:]])
        return combinations_arr

    # if mbdecomp is True, return array of every possible combination of index_list
    for n in range(1, len(index_list) + 1):
        # create array to hold all combinations of length n
        size_n_combinations = []
        for combination in itertools.combinations(index_list, n):
            size_n_combinations.append(combination)
        # append array of all combinations of size n to the main cominations array
        combinations_arr.append(size_n_combinations)
    return combinations_arr

def get_nmer_energies(molecule, config):
    """ 
    Computes the energy of a molecule and the MB decomposition,
    if requested.
    Input: A molecule that has no calculations done upon it yet, and a
           configuration file (settings.ini).
    Output: The calculated energy of the molecule. This can either have
            only one energy or many energies (from MB decomposition).
    """
    if config["MBdecomp"].getboolean("mbdecomp"):
        combinations = build_frag_indices(range(len(molecule.fragments)), True)
    else:
        combinations = build_frag_indices(range(len(molecule.fragments)), False)
    energy_str = ""
    output_str = ""
    for size_n_combinations in combinations:
        total_size_n_energy = []
        for combination in size_n_combinations:
            energy = calculator.calc_energy(molecule, combination, config)
            energy_str += str(energy) + " "
            output_str += "%.8f"%energy + " "
            molecule.energies[combination] = energy
            total_size_n_energy.append(energy)
        molecule.nmer_energies.append(total_size_n_energy)
    return output_str

# Get this printed to log file
# Loop through each k-mer, and for each k-mer do a k-body decomp
# Refine k-body energy function, we can discard some values returned from func
# This is called k-body decomp
def get_kbody_energies(mol):
    """
    Calculates the k-mer energies of the fragment
    subsets of a molecule. You can only run this function once you have
    calculated the many-body decomposition energies of the entire molecule.
    Input: A molecule with calculated many-body energies.
    Output: Many-body interaction energies of the subsets of molecules.
            Prints to terminal at the moment.
    """
    comb = build_frag_indices(range(len(mol.fragments)), len(mol.fragments))
    kbody_dict = {}
    kbody_diff = list(mol.mb_energies)
  
    # change comb to kbodys
    for kbody_subset in comb:
        for individual_kbody_indices in kbody_subset:
            kbody_indices = build_frag_indices(individual_kbody_indices, 
                len(individual_kbody_indices))
            kbody_energies = []
            for kbody_nfrags in kbody_indices:
                kbody_sub_arr = []
                for kbody_index in kbody_nfrags:
                    kbody_sub_arr.append(mol.energies[kbody_index])
                kbody_energies.append(kbody_sub_arr)
            kbody_input_index = kbody_indices[::-1]

            # Although mbdecomp() returns as an array of energies, we only
            # want the last one
            kbody_output_energy = mbdecomp(kbody_energies[::-1])[-1]

            # Likewise, we are only interested in one set of indices when
            # presenting our output
            kbody_output_index = kbody_input_index[-1][0]
            kbody_dict["K{}_{}B".format(kbody_output_index, 
                len(kbody_output_index))] = kbody_output_energy

            # Calculate k-body differences in comparison to our data
            which_body_energy = len(kbody_output_index) - 1
            kbody_diff[which_body_energy] -= kbody_output_energy
    return [kbody_dict, kbody_diff]

def mbdecomp(nmer_energies):
    """
    MB interaction energy decomposition
    Equation (5) from Gora et al., JCP 135 (2011) 224102
    Input: list of list of n-mer total energies,
           e.g. for a trimer [[1,2,3],[12,13,23],[123]]
    Output: list of nB interaction energies,
             e.g. for trimer [E1B,E2B,E3B]
    """
    
    # initialize array to hold return value
    nb_energies = []
    
    # i don't know what the rest of this does, Andy wrote it apparently
    esum = []
    nn = len(nmer_energies)
    for k in range(len(nmer_energies)):
        kk = k+1
        esum.append(sum(nmer_energies[k]))
        enb = 0
        for i in range(kk):
            ii = i+1
            a = nn - ii
            b = kk - ii
            fac = factorial(a)/(factorial(a-b)*factorial(b))
            if b % 2:
                fac = -1 * fac
            enb += fac * esum[i]
        nb_energies.append(enb)
    return nb_energies
