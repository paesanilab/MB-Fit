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

        self.variable_matcher = VariablePatternMatcher(variable_string)
        self.degree_matcher = NumberPatternMatcher(degree_string)

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

            if self.variable_matcher.match(variable.category) and self.degree_matcher.match(degree):
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

        self.variable_matcher = VariablePatternMatcher(variable_string)
        self.degree_matcher = NumberPatternMatcher(degree_string)

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

            # if this filter is  applicable to the current variable, add its degree to total_degree
            if self.variable_matcher.match(variable.category):
                total_degree += degree

        # now that we have summed the degrees of all filtered variables, check if that total degree is filtered
        # by this filter.
        if self.degree_matcher.match(total_degree):
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

        self.variable_matcher = VariablePatternMatcher(variable_string)
        self.fragment_num_matcher = NumberPatternMatcher(fragment_string)

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

            # if this filter is  applicable to the current variable, record the fragments it uses.
            if self.variable_matcher.match(variable.category):
                unique_fragments.add(variable.atom1_fragment)
                unique_fragments.add(variable.atom2_fragment)

        # count the number of fragments used by this monomial in all filtered variables
        num_fragments = len(unique_fragments)

        # now that we have found the number of fragments of all filtered variables, check if that number is filtered
        # by this filter.
        if self.fragment_num_matcher.match(num_fragments):
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

        self.variable_matcher = VariablePatternMatcher(variable_string)
        self.degree_matcher = NumberPatternMatcher(degree_string)
        self.term_matcher = NumberPatternMatcher(term_string)

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
        term = monomial.get_total_degree()

        # if this filter does not apply to the monomial's total degree, then it should not be filtered out
        if not self.term_matcher.match(term):
            return True

        # loop over each degree, variable pair in this monomial
        for degree, variable in zip(monomial, variables):
            
            if self.variable_matcher.match(variable.category) and self.degree_matcher.match(degree):
                return False

        # if the degree of one of this monomial's variables did not cause it to be filtered out, then it is a keeper!
        return True

class PatternMatcher(object):

    def match(self, string):
        raise NotImplementedError

class CompoundPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        patterns = pattern.split("/")

        self.pattern_matchers = []
        for pattern in patterns:
            self.pattern_matchers.append(self.get_sub_parser(pattern))

    def get_sub_parser(self, pattern):
        raise NotImplementedError

    def match(self, string):
        return any(pattern_matcher.match(string) for pattern_matcher in self.pattern_matchers)

class VariablePatternMatcher(CompoundPatternMatcher):

    def get_sub_parser(self, pattern):
        if pattern == "*":
            return WildCardPatternMatcher(pattern)
        else:
            var_type_pattern, atoms_pattern, level_pattern = pattern.split("-")[1:]

            if var_type_pattern == "*":
                var_type_matcher = WildCardPatternMatcher(var_type_pattern)
            else:
                var_type_matcher = StringPatternMatcher(var_type_pattern)

            if atoms_pattern == "*" or atoms_pattern == "*+*":
                atoms_matcher = WildCardPatternMatcher("*")
            elif atoms_pattern.split("+")[0] == "*":
                atoms_matcher = OneAtomMatcher(atoms_pattern.split("+")[1])
            elif atoms_pattern.split("+")[1] == "*":
                atoms_matcher = OneAtomMatcher(atoms_pattern.split("+")[0])
            else:
                atoms_matcher = BothAtomsMatcher(atoms_pattern)

            level_matcher = NumberPatternMatcher(level_pattern)

        return SingleVariablePatternMatcher(var_type_matcher, atoms_matcher, level_matcher)

class SingleVariablePatternMatcher(PatternMatcher):

    def __init__(self, var_type_matcher, atoms_matcher, level_matcher):
        self.var_type_matcher = var_type_matcher
        self.atoms_matcher = atoms_matcher
        self.level_matcher = level_matcher

    def match(self, string):
        var_type, atoms, level = string.split("-")[1:]

        return self.var_type_matcher.match(var_type) and self.atoms_matcher.match(atoms) and self.level_matcher.match(level)


class OneAtomMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.atom = pattern

    def match(self, string):
        atom1, atom2 = string.split("+")
        return self.atom == atom1 or self.atom == atom2

class BothAtomsMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.atom1, self.atom2 = pattern.split("+")

    def match(self, string):
        atom1, atom2 = string.split("+")
        return (self.atom1 == atom1 and self.atom2 == atom2) or (self.atom1 == atom2 and self.atom2 == atom1)

class StringPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.pattern = pattern

    def match(self, string):
        return self.pattern == string

class NumberPatternMatcher(CompoundPatternMatcher):

    def get_sub_parser(self, pattern):
        if pattern == "*":
            return WildCardPatternMatcher(pattern)
        elif pattern.endswith("-"):
            return LessThanOrEqualToPatternMatcher(pattern)
        elif pattern.endswith("+"):
            return GreaterThanOrEqualToPatternMatcher(pattern)
        elif "-" in pattern:
            return RangePatternMatcher(pattern)
        else:
            return EqualToPatternMatcher(pattern)

class LessThanOrEqualToPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.value = int(pattern[:-1])

    def match(self, string):
        return int(string) <= self.value

class GreaterThanOrEqualToPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.value = int(pattern[:-1])

    def match(self, string):
        return int(string) >= self.value

class EqualToPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.value = int(pattern)

    def match(self, string):
        return int(string) == self.value

class RangePatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        self.value1 = int(pattern[:pattern.index("-")])
        self.value2 = int(pattern[pattern.index("-") + 1:])

    def match(self, string):
        return int(string) >= self.value1 and int(string) <= self.value2

class WildCardPatternMatcher(PatternMatcher):

    def __init__(self, pattern):
        pass

    def match(self, string):
        return True
