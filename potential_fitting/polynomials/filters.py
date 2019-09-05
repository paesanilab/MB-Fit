# absolute module imports
from potential_fitting.exceptions import InvalidValueError, FilterBadSyntaxError

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

    input_args = args

    if len(args) == 0:
        raise FilterBadSyntaxError(input_args, 0, "nothing", "a valid filter name: 'not', 'degree', 'ind-degree', 'sum-degree', or 'num-fragments'")

    # keywords are not, and, or, degree
    if args[0] is '(':
        num_parens = 0
        for index, arg in enumerate(args):
            if arg is '(':
                num_parens += 1
            if arg is ')':
                num_parens -= 1
                if num_parens == 0:
                    break

        try:
            filter1 = parse_filter(*args[1:index])
        except FilterBadSyntaxError as e:
            raise FilterBadSyntaxError(input_args, 1 + e.index, e.saw, e.expected) from None

        args = args[index + 1:]
    else:
        # use the first argument to decide which filter to construct.
        if args[0] == "not":
            try:
                return NotFilter(parse_filter(*args[1:]))
            except FilterBadSyntaxError as e:
                raise FilterBadSyntaxError(input_args, 1 + e.index, e.saw, e.expected) from None
        if args[0] == "degree":
            try:
                filter1 = DegreeFilter(*args[1:4])
            except TypeError:
                raise FilterBadSyntaxError(input_args, len(args), "not enough arguments", "three arguments for 'degree' filter") from None
            index = 3
            args = args[index + 1:]
        elif args[0] == "ind-degree":
            try:
                filter1 = IndividualDegreeFilter(*args[1:3])
            except TypeError:
                raise FilterBadSyntaxError(input_args, len(args), "not enough arguments", "two arguments for 'ind-degree' filter") from None
            index = 2
            args = args[index + 1:]
        elif args[0] == "sum-degree":
            try:
                filter1 = SumDegreeFilter(*args[1:3])
            except TypeError:
                raise FilterBadSyntaxError(input_args, len(args), "not enough arguments", "two arguments for 'sum-degree' filter") from None
            index = 2
            args = args[index + 1:]
        elif args[0] == "num-fragments":
            try:
                filter1 = NumFragmentsFilter(*args[1:3])
            except TypeError:
                raise FilterBadSyntaxError(input_args, len(args), "not enough arguments", "two arguments for 'num-fragments' filter") from None
            index = 2
            args = args[index + 1:]
        else:
            raise FilterBadSyntaxError(input_args, 0, args[0], "a valid filter name: 'not', 'degree', 'ind-degree', 'sum-degree', or 'num-fragments'")


    if len(args) == 0:
        return filter1

    if args[0] != "and" and args[0] != "or":
        raise FilterBadSyntaxError(input_args, index + 1, args[0], "either 'and' or 'or'")

    conjunction = args[0]

    args = args[1:]

    if len(args) == 0:
        raise FilterBadSyntaxError(input_args, index + 2, "no statement", "a statement after conjunction '{}'".format(conjunction))

    try:
        filter2 = parse_filter(*args)
    except FilterBadSyntaxError as e:
        raise FilterBadSyntaxError(input_args, index + 2 + e.index, e.saw, e.expected) from None

    if conjunction == "and":
        return AndFilter(filter1, filter2)
    if conjunction == "or":
        return OrFilter(filter1, filter2)

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

class AndFilter(Filter):
    """
    Combines two filters, so this filter will filter out any terms that would be filtered out by both other filters.
    """

    def __init__(self, filter1, filter2):
        """
        Creates a new AndFilter from two filters to combine.

        This filter will filter out terms that are filtered out by both filter1 and filter2.

        Args:
            filter1         - The first filter to combine.
            filter2         - The second filter to combine.
        """
        self.filter1 = filter1
        self.filter2 = filter2

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # if this monomial is kept by either filter1 or filter2, then it is not filtered out by both of them,
        # so keep this monomial.
        return self.filter1.keep(monomial, variables) or self.filter2.keep(monomial, variables)

class OrFilter(Filter):
    """
    Combines two filters, so this filter will filter out any terms that would be filtered out by either other filter.
    """

    def __init__(self, filter1, filter2):
        """
        Creates a new AndFilter from two filters to combine.

        This filter will filter out terms that are filtered out by either filter1 or filter2.

        Args:
            filter1         - The first filter to combine.
            filter2         - The second filter to combine.
        """
        self.filter1 = filter1
        self.filter2 = filter2

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # if this monomial is kept by both filter1 and filter2, then it is not filtered out by either of them,
        # so keep this monomial.
        return self.filter1.keep(monomial, variables) and self.filter2.keep(monomial, variables)

class IndividualDegreeFilter(Filter):
    """
    Filters out monomials based on their individual degree in one or more particular variables.
    """

    def __init__(self, variable_string, degree_string):
        """
        Creates a new IndividualDegreeFilter from the given parameters.

        This filter will filter OUT any monomials with degree as specified by degree_string in one of the variables as
        specified by variable_string.

        Each variable that matches variable_string is treated independantly, meaning they are not summed.

        Args:
            variable_string - '/' delimited list of variables to apply this filter to. Each item seperated by a '/'
                    should be in the following format:
                        * TYPE-ATOMS
                    where TYPE is one of:
                        * x         -- Only affects intermolecular variables.
                        * x-intra   -- Only affects intramolecular variables.
                        * x-*       -- Affects both intermolecular and intramolecular variables.
                    and ATOMS is one of:
                        * A+B        -- Affects variables describing the distance between A and B.
                        * *+A or A+*  -- Affects variables describing the distance between A and any other atom.
                        * *+*        -- Affects variables regardless of the atom types involved.
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

        Returns:
            A new IndividualDegreeFilter.

        """

        self.variables = variable_string.split("/")
        self.degrees = degree_string.split("/")

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

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

                    var_type, atoms_string = variable_string.split("-")[1:]
                    atom1, atom2 = atoms_string.split("+")

                    # if the variable string ends in 2 wildcards, then this filter is eligable to apply to all
                    # variables regardless of atom type
                    if atom1 == "*" and atom2 == "*":
                        fits_atoms = True

                    # if the second to last character or the last character of the variable string is a wildcard (but
                    # not both), then this filter is eligable to apply to all variables where one of the atoms is the
                    # other atom specified by the variable string
                    elif atom1 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom2
                                      or variable.category.split("-")[-1].split("+")[1] == atom2)
                    elif atom2 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom1
                                      or variable.category.split("-")[-1].split("+")[1] == atom1)

                    # otherwise, this variable string simply specifies an exact pair of atoms to be effected by this
                    # filter
                    else:
                        fits_atoms = variable.category.split("-")[-1].split("+")[0] == atom1 and \
                                     variable.category.split("-")[-1].split("+")[1] == atom2

                    # check if the middle bit of this filter is an "*"
                    if variable_string.split("-")[1] == "*":
                        fits_type = True

                    # otherwise, this variable string specifies an exact inter/intra type (x- or x-intra-)
                    else:
                        fits_type = variable_string.split("-")[1] == variable.category.split("-")[1]

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
                    if (degree >= int(degree_string[:degree_string.index("-")]) and degree <=
                            int(degree_string[degree_string.index("-") + 1:])):
                        return False

                # otherwise, it is an exact value specifier
                else:
                    if degree == int(degree_string):
                        return False

        # if the degree of one of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True

class SumDegreeFilter(Filter):
    """
    Filters out monomials based on their total degree in one or more particular variables.
    """

    def __init__(self, variable_string, degree_string):
        """
        Creates a new SumDegreeFilter from the given parameters.

        This filter will filter OUT any monomials with total degree as specified by degree_string in all of the variables as
        specified by variable_string.

        The degree of all variables matching variable_string are summed, and the compared to degree_string.

        Args:
            variable_string - '/' delimited list of variables to apply this filter to. Each item seperated by a '/'
                    should be in the following format:
                        * TYPE-ATOMS
                    where TYPE is one of:
                        * x         -- Only affects intermolecular variables.
                        * x-intra   -- Only affects intramolecular variables.
                        * x-*       -- Affects both intermolecular and intramolecular variables.
                    and ATOMS is one of:
                        * A+B        -- Affects variables describing the distance between A and B.
                        * *+A or A+*  -- Affects variables describing the distance between A and any other atom.
                        * *+*        -- Affects variables regardless of the atom types involved.
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

        Returns:
            A new SumDegreeFilter.

        """

        self.variables = variable_string.split("/")
        self.degrees = degree_string.split("/")

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # keep track of the total degree in all applicable variables
        total_degree = 0

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

                    var_type, atoms_string = variable_string.split("-")[1:]
                    atom1, atom2 = atoms_string.split("+")

                    # if the variable string ends in 2 wildcards, then this filter is eligable to apply to all
                    # variables regardless of atom type
                    if atom1 == "*" and atom2 == "*":
                        fits_atoms = True

                    # if the second to last character or the last character of the variable string is a wildcard (but
                    # not both), then this filter is eligable to apply to all variables where one of the atoms is the
                    # other atom specified by the variable string
                    elif atom1 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom2
                                      or variable.category.split("-")[-1].split("+")[1] == atom2)
                    elif atom2 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom1
                                      or variable.category.split("-")[-1].split("+")[1] == atom1)

                    # otherwise, this variable string simply specifies an exact pair of atoms to be effected by this
                    # filter
                    else:
                        fits_atoms = variable.category.split("-")[-1].split("+")[0] == atom1 and \
                                     variable.category.split("-")[-1].split("+")[1] == atom2

                    # check if the middle bit of this filter is an "*"
                    if variable_string.split("-")[1] == "*":
                        fits_type = True

                    # otherwise, this variable string specifies an exact inter/intra type (x- or x-intra-)
                    else:
                        fits_type = variable_string.split("-")[1] == variable.category.split("-")[1]

                    # if this variable string fits both the atoms and the inter/intra type, then it is applicable to
                    # this variable!
                    if fits_type and fits_atoms:
                        applicable_variable = True
                        break

                # otherwise, this string simply specifies an exact variable
                elif variable_string == variable.category:
                    applicable_variable = True
                    break

            # if this filter is  applicable to the current variable, add its degree to total_degree
            if applicable_variable:
                total_degree += degree

        # now that we have summed the degrees of all filtered variables, check if that total degree is filtered
        # by this filter.

        # loop over each degree affected by this filter
        for degree_string in self.degrees:

            # if this degree string is the wildcard, then it affects all variables regardless of their degree
            if degree_string == "*":
                return False

            # if this degree string ends in a "-", then it is a less than or equal to specifier
            elif degree_string.endswith("-"):
                if total_degree <= int(degree_string[:-1]):
                    return False

            # if this degree string ends in a "+", then it is a greater than or equal to specifier
            elif degree_string.endswith("+"):
                if total_degree >= int(degree_string[:-1]):
                    return False

            # if this degree contains a "-" (but does not end in "-") then it is a range specifier
            elif "-" in degree_string:
                if (total_degree >= int(degree_string[:degree_string.index("-")]) and degree <=
                        int(degree_string[degree_string.index("-") + 1:])):
                    return False

            # otherwise, it is an exact value specifier
            else:
                if total_degree == int(degree_string):
                    return False

        # if the total_degree of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True

class NumFragmentsFilter(Filter):
    """
    Filters out monomials based on the number of fragments in a monomial.
    """

    def __init__(self, variable_string, fragment_string):
        """
        Creates a new NumFragmentsFilter from the given parameters.

        This filter will filter OUT any monomials with total fragments as specified by fragment_string in all of the variables as
        specified by variable_string.

        A variable is considered to include a fragment if one if its atoms are in that fragment.

        The number of fragments of all variables matching variable_string are summed, and the compared to fragment_string.

        Args:
            variable_string - '/' delimited list of variables to apply this filter to. Each item seperated by a '/'
                    should be in the following format:
                        * TYPE-ATOMS
                    where TYPE is one of:
                        * x         -- Only affects intermolecular variables.
                        * x-intra   -- Only affects intramolecular variables.
                        * x-*       -- Affects both intermolecular and intramolecular variables.
                    and ATOMS is one of:
                        * A+B        -- Affects variables describing the distance between A and B.
                        * *+A or A+*  -- Affects variables describing the distance between A and any other atom.
                        * *+*        -- Affects variables regardless of the atom types involved.
                    A and B can be substituted for any other atom type.
            fragment_string   - '/' delimited list of fragment numbers to apply this filter to. Each item seperated by a '/'
                    should be in one of the following formats:
                        * y-        -- Affects monomials with total fragments equal to or less than y in one or more of the
                                specified variables.
                        * y+        -- Affects monomials with total fragments equal to or greater than y in one or more of the
                                specified variables.
                        * y-z       -- Affects monomials with total fragments in the range [y, z] (inclusive) in one or more of
                                the specified variables.
                        * y         -- Affects monomials with total fragments y in one or more of the specified variables.
                        * *         -- Affects monomials regardless of their degrees.
                    y and z can be any integer >= 0.

        Returns:
            A new NumFragmentsFilter.

        """

        self.variables = variable_string.split("/")
        self.fragment_nums = fragment_string.split("/")

    def keep(self, monomial, variables):
        """
        Tells whether the input monomial formed by the input variables is not filtered out by this filter.

        Args:
            monomial        - The monomial to filter, specified as list of degrees of length len(variables).
            variables       - List of variables in this monomial, should be same length as monomial.

        Returns:
            False if this Filter filters out this monomial, True otherwise.
        """

        # keep track of all fragments used in this monomial
        unique_fragments = set()

        # loop over each degree, variable pair in this monomial
        for degree, variable in zip(monomial, variables):

            # if this variable is not used in the monomial, do not add its fragments to unique_fragments.
            if degree == 0:
                continue

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

                    var_type, atoms_string = variable_string.split("-")[1:]
                    atom1, atom2 = atoms_string.split("+")

                    # if the variable string ends in 2 wildcards, then this filter is eligable to apply to all
                    # variables regardless of atom type
                    if atom1 == "*" and atom2 == "*":
                        fits_atoms = True

                    # if the second to last character or the last character of the variable string is a wildcard (but
                    # not both), then this filter is eligable to apply to all variables where one of the atoms is the
                    # other atom specified by the variable string
                    elif atom1 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom2
                                      or variable.category.split("-")[-1].split("+")[1] == atom2)
                    elif atom2 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom1
                                      or variable.category.split("-")[-1].split("+")[1] == atom1)

                    # otherwise, this variable string simply specifies an exact pair of atoms to be effected by this
                    # filter
                    else:
                        fits_atoms = variable.category.split("-")[-1].split("+")[0] == atom1 and \
                                     variable.category.split("-")[-1].split("+")[1] == atom2

                    # check if the middle bit of this filter is an "*"
                    if variable_string.split("-")[1] == "*":
                        fits_type = True

                    # otherwise, this variable string specifies an exact inter/intra type (x- or x-intra-)
                    else:
                        fits_type = variable_string.split("-")[1] == variable.category.split("-")[1]

                    # if this variable string fits both the atoms and the inter/intra type, then it is applicable to
                    # this variable!
                    if fits_type and fits_atoms:
                        applicable_variable = True
                        break

                # otherwise, this string simply specifies an exact variable
                elif variable_string == variable.category:
                    applicable_variable = True
                    break

            # if this filter is  applicable to the current variable, record the fragments it uses.
            if applicable_variable:
                unique_fragments.add(variable.atom1_fragment)
                unique_fragments.add(variable.atom2_fragment)

        # count the number of fragments used by this monomial in all filtered variables
        num_fragments = len(unique_fragments)

        # now that we have found the number of fragments of all filtered variables, check if that number is filtered
        # by this filter.

        # loop over each fragment number affected by this filter
        for fragment_string in self.fragment_nums:

            # if this fragment_num string is the wildcard, then it affects all variables regardless of their fragment number
            if fragment_string == "*":
                return False

            # if this fragment_num string ends in a "-", then it is a less than or equal to specifier
            elif fragment_string.endswith("-"):
                if num_fragments <= int(fragment_string[:-1]):
                    return False

            # if this fragment_num string ends in a "+", then it is a greater than or equal to specifier
            elif fragment_string.endswith("+"):
                if num_fragments >= int(fragment_string[:-1]):
                    return False

            # if this fragment_num contains a "-" (but does not end in "-") then it is a range specifier
            elif "-" in fragment_string:
                if (num_fragments >= int(fragment_string[:fragment_string.index("-")]) and degree <=
                        int(fragment_string[fragment_string.index("-") + 1:])):
                    return False

            # otherwise, it is an exact value specifier
            else:
                if num_fragments == int(fragment_string):
                    return False

        # if the fragment number of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True

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
                        * A+B        -- Affects variables describing the distance between A and B.
                        * *+A or A+*  -- Affects variables describing the distance between A and any other atom.
                        * *+*        -- Affects variables regardless of the atom types involved.
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

                    var_type, atoms_string = variable_string.split("-")[1:]
                    atom1, atom2 = atoms_string.split("+")

                    # if the variable string ends in 2 wildcards, then this filter is eligable to apply to all
                    # variables regardless of atom type
                    if atom1 == "*" and atom2 == "*":
                        fits_atoms = True

                    # if the second to last character or the last character of the variable string is a wildcard (but
                    # not both), then this filter is eligable to apply to all variables where one of the atoms is the
                    # other atom specified by the variable string
                    elif atom1 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom2
                                or variable.category.split("-")[-1].split("+")[1] == atom2)
                    elif atom2 == "*":
                        fits_atoms = (variable.category.split("-")[-1].split("+")[0] == atom1
                                or variable.category.split("-")[-1].split("+")[1] == atom1)

                    # otherwise, this variable string simply specifies an exact pair of atoms to be effected by this
                    # filter
                    else:
                        fits_atoms = variable.category.split("-")[-1].split("+")[0] == atom1 and variable.category.split("-")[-1].split("+")[1] == atom2

                    # check if the middle bit of this filter is an "*"
                    if variable_string.split("-")[1] == "*":
                        fits_type = True

                    # otherwise, this variable string specifies an exact inter/intra type (x- or x-intra-)
                    else:
                        fits_type = variable_string.split("-")[1] == variable.category.split("-")[1]

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
                    if (degree >= int(degree_string[:degree_string.index("-")]) and degree <=
                            int(degree_string[degree_string.index("-") + 1:])):
                        return False

                # otherwise, it is an exact value specifier
                else:
                    if degree == int(degree_string):
                        return False

        # if the degree of one of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True
