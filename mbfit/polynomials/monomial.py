import hashlib

class Monomial(object):
    """
    Class that stores all the information associated with a single Monomial: one term of a polynomial.
    """

    def __init__(self, degrees):
        """
        Constructs a new Monomial.

        Args:
            degrees     - List of the degrees of each variable in this monomial.

        Returns:
            A new Monomial.
        """

        self.degrees = degrees

    def permute(self, variable_permutations):
        """
        Gives all permutations of this Monomial.

        Args:
            variable_permutations - Variable permutations as generated by get_variable_permutations().

        Generates:
            List of all Monomials created by permuting the this Monomial by all the variable_permutations.
        """

        # there will be 1 permutation for each variable permutation.
        for variable_permutation in variable_permutations:

            # initialize the permutation to all zeros
            monomial_permutation = [0 for i in self.degrees]

            # loop over each variable index and degree
            for index, degree in zip(range(len(self.degrees)), self.degrees):

                # because we initialized the monomial_permutation to all 0s, if degree is 0, we don't have to do anything
                if degree == 0:
                    continue

                # if degree is not zero, then lookup the new index of this variable in the permuted monomial
                new_index = variable_permutation.index(index)

                # now set the value of this variable in the permuted monomial
                monomial_permutation[new_index] = self.degrees[index]

            yield Monomial(monomial_permutation)

    def get_standard_permutations(self, variable_permutations):
        return sorted(self.permute(variable_permutations), key=lambda x: x.get_degrees())[-1]

    def get_total_degree(self):
        """
        Gets the total degree of this Monomial by summing the degree of each variable.

        Args:
            None.

        Returns:
            The total degree of this Monomial.
        """
        return sum(self.degrees)

    def get_derivative(self, variable_index):
        """
        Gets the derivative of this Monomial with respect to the variable at the given index.

        Args:
            variable_index      - The index of the variable to take the derivative with respect to.

        Returns:
            (constant, monomial)
            constant - a constant that the derivative must be multiplied by.
            monomial - A new Monomial, such that the derivative of this Monomial with respect to the variable at the
                    given index is constant*monomial.
        """

        # allow the user to specify distance from end of the Monomial by giving a negative variable_index.
        if variable_index < 0:
            variable_index = len(self.degrees) + variable_index

        return self.degrees[variable_index], Monomial([d - 1 if i == variable_index else d
                                                       for i, d
                                                       in enumerate(self.degrees)])

    def get_degrees(self):
        return self.degrees

    def __eq__(self, other):
        return self.degrees == other.degrees

    def __hash__(self):
        return int(hashlib.sha1(":".join([str(degree) for degree in self.degrees]).encode()).hexdigest(), 16)

    def __getitem__(self, item):
        return self.degrees[item]
