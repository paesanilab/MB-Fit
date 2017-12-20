import itertools

def make_combs(size):
    """ Make an array of all combinations possible.
        Note that the first value of the list helps
        with total energy.
        Input: The largest size of the combination.
        Output: An array of combinations of ints.
    """
    comb_set = range(size)
    comb_arr = []

    while size > 0:
        for comb in itertools.combinations(comb_set, size):
            comb_arr.append(comb)
        size -= 1
    return comb_arr
