# external package imports
import itertools
from math import factorial

# absolute module imports
from potential_fitting.utils import SettingsReader

# local module imports
from .calculator_utils import get_calculator
from .model import Model

# ****** I'm not gonna clean this up because its on the top of the rework/remove chopping block ******
# ****** currently only used by driver.py ******

def build_frag_indices(index_list, mbdecomp):
    """
    Returns a 3d list of indicies to be passed into calc_energy.
    Its hard to explain what this does, so I'll just give examples:    

    if mbdecomp is true, returns an list of every possible combinations
    in the list, sorted into list by size.
        Example: index_list = [0, 1, 2]
                return value = [
                    [[0], [1], [2]],          # list of size 1 combinations
                    [[0, 1], [0, 2], [1, 2]], # list of size 2 combinations
                    [[0, 1, 2]]               # list of size 3 combinations
                ]

    if mbdecomp is false:
        Example: index_list = [0, 1, 2]
                return value = [[[0, 1, 2]]]  # list with only size 3 combination
    """
    # used to hold list of combinations of index_list
    combinations_arr = []
    
    # if mbdecomp is False, just return a copy of index_list, inside a two additional layers of list
    if mbdecomp == False:
        combinations_arr.append([index_list[:]])
        return combinations_arr

    # if mbdecomp is True, return list of every possible combination of index_list
    for n in range(1, len(index_list) + 1):
        # create list to hold all combinations of length n
        size_n_combinations = []
        for combination in itertools.combinations(index_list, n):
            size_n_combinations.append(combination)
        # append list of all combinations of size n to the main cominations list
        combinations_arr.append(size_n_combinations)
    return combinations_arr

def get_nmer_energies(settings_path, molecule):

    """ 
    Computes the energy of a molecule and the MB decomposition,
    if requested.
    Input: A molecule that has no calculations done upon it yet, and a
           configuration file (settings.ini).
    Output: The calculated energy of the molecule. This can either have
            only one energy or many energies (from MB decomposition).
    """
    settings = SettingsReader(settings_path)
    if settings.getboolean("MBdecomp", "mbdecomp"):
        # if multibody decomposition is requested, get a list of all possible combinations of fragment indices sorted into sub lists by size
        combinations = build_frag_indices(range(len(molecule.fragments)), True)
    else:
        # if not, the list will only contain 1 sublist with 1 combination
        combinations = build_frag_indices(range(len(molecule.fragments)), False)
    
    # initialize string to build return value
    output_str = ""

    # for every set of size n combinations...
    for size_n_combinations in combinations:
        # list of all energies of fragment combinations of size n
        total_size_n_energy = []
        
        # for each combination of size n...
        for combination in size_n_combinations:
            # calculate the energy of this set of fragments
            calc = get_calculator(settings_path, logging=True)
            energy, log_path = calc.calculate_energy(molecule, Model(settings.get("model", "method"),
                                                                     settings.get("model", "basis"),
                                                                     settings.getboolean("model", "cp")), combination)
            # build output_str with 9
            output_str += "%.8f"%energy + " "
            molecule.energies[combination] = energy
            total_size_n_energy.append(energy)
        molecule.nmer_energies.append(total_size_n_energy)
    return output_str

# Get this printed to log file
# Loop through each k-mer, and for each k-mer do a k-body decomp
# Refine k-body energy function, we can discard some values returned from func
# This is called k-body decomp
def get_kbody_energies(molecule):
    """
    Calculates the k-mer energies of the fragment
    subsets of a molecule. You can only run this function once you have
    calculated the many-body decomposition energies of the entire molecule.
    Input: A molecule with calculated many-body energies.
    Output: Many-body interaction energies of the subsets of molecules.
            Prints to terminal at the moment.
    """

    # get list of every possibe combinations of the indicies of the fragments,
    # sorted by size. ie: [[1,2,3],[12,23,13],[123]]
    nbody_combinations = build_frag_indices(range(len(molecule.fragments)), True)
    # 
    kbody_dict = {}
    # 
    kbody_diff = list(molecule.mb_energies)
  
    
    # for each list of size n combinations of fragments
    for size_n_combinations in nbody_combinations:
        # for each combination of size n
        for size_n_combination in size_n_combinations:
            # create list of all possible combinations of that combination
            kbody_combinations = build_frag_indices(size_n_combination, True)
            # list of kbody energies, sorted into sublists by size. ie: [[1, 2, 3], [12,23,13], [123]]
            kbody_energies = []
            # for each list of size k combinations of fragments
            for size_k_combinations in kbody_combinations:
                # list of energies for size k bodies
                size_k_energies = []
                # for each combination of size k
                for size_k_combination in size_k_combinations:
                    # get the energy of that combination
                    size_k_energies.append(molecule.energies[size_k_combination])
                # add the list of energies for all combinations of size k to the total list of energies
                kbody_energies.append(size_k_energies)

            kbody_input_index = kbody_combinations

            # Although mbdecomp() returns as an list of energies, we only
            # want the last one
            kbody_output_energy = mbdecomp(kbody_energies)[-1]

            # Likewise, we are only interested in one set of indices when
            # presenting our output
            kbody_output_index = kbody_input_index[-1][0]
            kbody_dict["V{}_{}B".format(kbody_output_index, 
                len(kbody_output_index))] = kbody_output_energy

            # Calculate k-body differences in comparison to our data
            which_body_energy = len(kbody_output_index) - 1
            kbody_diff[which_body_energy] -= kbody_output_energy
    return [kbody_dict, kbody_diff]


"""
Magic
"""
def mbdecomp(nmer_energies):
    """
    MB interaction energy decomposition
    Equation (5) from Gora et al., JCP 135 (2011) 224102
    Input: list of list of n-mer total energies,
           e.g. for a trimer [[1,2,3],[12,13,23],[123]]
    Output: list of nB interaction energies,
             e.g. for trimer [E1B,E2B,E3B]
    """
    
    # initialize list to hold return value
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
