# absolute module imports
from potential_fitting.exceptions import InvalidValueError

def parse_filter(*args):
    """
    Takes in the arguments for a single filter and returns the filter that they represent.

    Valid filters are Not and Degree, see descriptions of their constructors to see their arguments.

    Args:
        args                - The arguments that make up the filter. The first item should be the name of the filter
                ("not" or "degree") the rest should be the arguments for the constructor of that filter.

    Returns:
        A new filter object.
    """

    # use the first argument to decide which filter to construct.
    if args[0] == "not":
        return NotFilter(parse_filter(*args[1:]))
    if args[0] == "degree":
        return DegreeFilter(*args[1:])

    # if the filter type is not a valid filter, raise an error.
    raise InvalidValueError("filter type", args[0], "must be one of 'not', 'degree'")

class Filter(object):
    """
    Abstract Class for all filters to extend.
    """

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        raise NotImplementedError

class NotFilter(Filter):
    """
    Inverts another filter, so this filter will filter out any terms that would NOT be filtered out by the other
    filter.
    """

    def __init__(self, not_filter):
        """
        Creates a new NotFilter from a filter to invert.

        Args:
            not_filter      - This filter is the filter to be inverted by this one.
        """
        self.not_filter = not_filter

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # simply check if this monomial would not be filtered out by the inverted filter
        return not self.not_filter.keep(monomial, variables)


class DegreeFilter(Filter):
    """
    Filters out monomials based on their degree in particular variables and their overal degree.
    """

    def __init__(self, variable_string, degree_string, term_string):
        """
        Creates a new DegreeFilter from the given parameters.

        This filter will filter OUT any monomials with degree as specified by degree_string in one of the variables as
        specified by variable_string with TOTAL degree as specified by term_string.

        Args:
            variable_string - '/' delimited list of variables to apply this filter to. Each item seperated by a '/'
                    should be in the following format:
                        * TYPE-ATOMS
                    where TYPE is one of:
                        * x         -- Only affects intermolecular variables.
                        * x-intra   -- Only affects intramolecular variables.
                        * x-*       -- Affects both intermolecular and intramolecular variables.
                    and ATOMS is one of:
                        * AB        -- Affects variables describing the distance between A and B.
                        * *A or A*  -- Affects variables describing the distance between A and any other atom.
                        * **        -- Affects variables regardless of the atom types involved.
                    A and B can be substituted for any other atom type.
            degree_string   - '/' delimited list of degrees to apply this filter to. Each item seperated by a '/'
                    should be in one of the following formats:
                        * y-        -- Affects monomials with degree equal to or less than y in one or more of the
                                specified variables.
                        * y+        -- Affects monomials with degree equal to or greater than y in one or more of the
                                specified variables.
                        * y-z       -- Affects monomials with degree in the range [y, z] (inclusive) in one or more of
                                the specified variables.
                        * y         -- Affects monomials with degree y in one or more of the specified variables.
                        * *         -- Affects monomials regardless of their degrees.
                    y and z can be any integer >= 0.
            term_string - '/' delimited list of terms to apply this filter to. Each item seperated by a '/' should be
                    in one of the following formats:
                        * y-        -- Affects monomials with TOTAL degree equal to or less than y.
                        * y+        -- Affects monomials with TOTAL degree equal to or greater than y.
                        * y-z       -- Affects monomials with TOTAL degree in the range [y, z] (inclusive).
                        * y         -- Affects monomials with TOTAL degree y.
                        * *         -- Affects monomials regardless of their TOTAL degree.
                    y and z can be any integer >= 0.

        Returns:
            A new DegreeFilter.

        """

        self.variables = variable_string.split("/")
        self.degrees = degree_string.split("/")
        self.terms = term_string.split("/")

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # calculate the total degree of this monomial
        term = sum(monomial)

        # used to track whether this monomial's TOTAl degree matches one of this filter's applicable terms
        applicable_term = False

        # loop over each term string this filter effects
        for term_string in self.terms:

            # if the term string is the wildcard, this filter is applicable to all monomials
            if term_string == "*":
                applicable_term = True
                break

            # if the term string ends in "-" then it is a less than or equal to specifier
            elif term_string.endswith("-"):
                if term <= int(term_string[:-1]):
                    applicable_term = True
                    break

            # if the term string ends in "+" then it is a greater than or equal to specifier
            elif term_string.endswith("+"):
                if term >= int(term_string[:-1]):
                    applicable_term = True
                    break

            # if the term has a "-" (but did not end in "-") then it is a range specifier
            elif "-" in term_string:
                if term >= int(term_string[:term_string.index("-")]) and term <= int(term_string[term_string.index("-")
                        + 1:]):
                    applicable_term = True
                    break

            # otherwise, it is a specific number specifier
            else:
                if term == int(term_string):
                    applicable_term = True
                    break

        # if this filter does not apply to the monomial's total degree, then it should not be filtered out
        if not applicable_term:
            return True

        # loop over each degree, variable pair in this monomial
        for degree, variable in zip(monomial, variables):

            # used to track whether this variable matches one of this filter's applicable variables
            applicable_variable = False

            # loop over each variable this filter effects
            for variable_string in self.variables:

                # if the variable string is the wildcard, this filter is applicable to all variables
                if variable_string == "*":
                    applicable_variable = True
                    break

                # if the variable string has an *, we must perform detailed analysis to figure out what it does
                if "*" in variable_string:

                    # tracks if the atom portion of this variable string matches the current variable's
                    fits_atoms = False
                    # tracks if the type portion of this variable string matches the current variable's (inter or
                    # intra)
                    fits_type = False

                    # if the variable string ends in 2 wildcards, then this filter is eligable to apply to all
                    # variables regardless of atom type
                    if variable_string[-2] == "*" and variable_string[-1] == "*":
                        fits_atoms = True

                    # if the second to last character or the last character of the variable string is a wildcard (but
                    # not both), then this filter is eligable to apply to all variables where one of the atoms is the
                    # other atom specified by the variable string
                    elif variable_string[-2] == "*":
                        fits_atoms = (variable.category[-1] == variable_string[-1]
                                or variable.category[-2] ==variable_string[-1])
                    elif variable_string[-1] == "*":
                        fits_atoms = (variable.category[-1] == variable_string[-2]
                                or variable.category[-2] == variable_string[-2])

                    # otherwise, this variable string simply specifies an exact pair of atoms to be effected by this
                    # filter
                    else:
                        fits_atoms = variable.category[-2:] == variable_string[-2:]

                    # if the second character of the variable string is an * but not the 3rd (as in x-*-AA but NOT
                    # x-**), then this filter is eligable to apply to all variables regardless of inter/intra type
                    if variable_string[2] == "*" and variable_string[3] != "*":
                        fits_type = True

                    # otherwise, this variable string specifies an exact inter/intra type (x- or x-intra-)
                    else:
                        fits_type = variable_string[:-2] == variable.category[:-2]

                    # if this variable string fits both the atoms and the inter/intra type, then it is applicable to
                    # this variable! 
                    if fits_type and fits_atoms:
                        applicable_variable = True
                        break 

                # otherwise, this string simply specifies an exact variable
                elif variable_string == variable.category:
                    applicable_variable = True
                    break
            
            # if this filter is not applicable to the current variable, go to the next one
            if not applicable_variable:
                continue

            # loop over each degree affected by this filter
            for degree_string in self.degrees:

                # if this degree string is the wildcard, then it affects all variables regardless of their degree
                if degree_string == "*":
                    return False

                # if this degree string ends in a "-", then it is a less than or equal to specifier
                elif degree_string.endswith("-"):
                    if degree <= int(degree_string[:-1]):
                        return False

                # if this degree string ends in a "+", then it is a greater than or equal to specifier
                elif degree_string.endswith("+"):
                    if degree >= int(degree_string[:-1]):
                        return False

                # if this degree contains a "-" (but does not end in "-") then it is a range specifier
                elif "-" in degree_string:
                    if degree >= int(degree_string[:degree_string.index("-")]) and degree <=
                            int(degree_string[degree_string.index("-") + 1:]):
                        return False

                # otherwise, it is an exact value specifier
                else:
                    if degree == int(degree_string):
                        return False

        # if the degree of one of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True


