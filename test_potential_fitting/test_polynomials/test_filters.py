import unittest, random

from potential_fitting.polynomials import filters, Variable, Monomial
from potential_fitting.exceptions import FilterBadSyntaxError

class TestFilters(unittest.TestCase):

    def test_parse_no_arguments(self):

        # Test for when you give parse_filter() no arguments

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter()

    def test_degree_wrong_num_args(self):

        # Tests for when you give 'degree' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', 'arg1')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', 'arg1', 'arg2')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', 'arg1', 'arg2', 'arg3', 'arg4')

    def test_ind_degree_wrong_num_args(self):

        # Tests for when you give 'ind-degre' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', 'arg1')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', 'arg1', 'arg2', 'arg3')

    def test_sum_degree_wrong_num_args(self):

        # Tests for when you give 'sum-degree' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree', 'arg1')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree', 'arg1', 'arg2', 'arg3')

    def test_num_fragments_wrong_num_args(self):

        # Tests for when you give 'num-fragments' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments', 'arg1')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments', 'arg1', 'arg2', 'arg3')

    def test_unrecognized_filter_name(self):

        # Test for unrecognized filter name

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('some-bad-filter')

    def test_dangling_conjunction(self):

        # Tests for dangling conjunction

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', 'arg1', 'arg2', 'and')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', 'arg1', 'arg2', 'or')

    def test_parse_dangling_not(self):

        # Tests for danging not
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not')

    def test_parse_bad_syntax_in_parens(self):

        # Test for when bad syntax appears within parenthesis

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'ind-degree', 'arg1', 'arg2', 'and', ')')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'ind-degree', 'arg1', ')')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'some-bad-filter', ')')

    def test_parse_bad_syntax_in_not(self):

        # Test for when bad syntax appears within a not

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'ind-degree', 'arg1')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'not')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'some-bad-filter')

    def test_parse_bad_syntax_after_conjunction(self):

        # Test for when bad syntax appears after a conjunction

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'and', 'sum-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'and', 'some-bad-filter')

    def test_filter(self):
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]
        with self.assertRaises(NotImplementedError):
            filters.Filter().keep(Monomial([0, 1, 1]), variables)

    def test_parse_individual_degree_filter(self):
        filter = filters.parse_filter('ind-degree', "*", "2+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

    def test_parse_sum_degree_filter(self):
        filter = filters.parse_filter("sum-degree", "*", "2+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

    def test_parse_degree_filter(self):
        filter = filters.parse_filter("degree", "*", "2+", "3")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([15, 2, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))

    def test_parse_num_fragments_filter(self):
        filter = filters.parse_filter("num-fragments", "*", "1-/3+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'b', 'A3', 'b', 'x-intra-A+A'),
                     Variable('A2', 'c', 'A3', 'd', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([3, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 3]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))

    def test_parse_not_filter(self):
        filter = filters.IndividualDegreeFilter("*", "2+")
        not_filter = filters.parse_filter("not", "ind-degree", "*", "2+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        filter = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        not_filter = filters.parse_filter("not", "num-fragments", "x-intra-A+B", "1-")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

    def test_parse_and_filter(self):
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        and_filter = filters.parse_filter("ind-degree", "*", "2+", "and", "num-fragments", "x-intra-A+B", "1-")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) or filter2.keep(monomial, variables),
                            and_filter.keep(monomial, variables))

    def test_parse_or_filter(self):
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        or_filter = filters.parse_filter("ind-degree", "*", "2+", "or", "num-fragments", "x-intra-A+B", "1-")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) and filter2.keep(monomial, variables),
                             or_filter.keep(monomial, variables))

    def test_individual_degree_filter(self):
        filter = filters.IndividualDegreeFilter("*", "2+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+B", "1/3+")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'a', 'B2', 'a', 'x-intra-A+B'),
                     Variable('B1', 'a', 'B2', 'a', 'x-intra-B+B')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 5]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 2, 3]), variables))
        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 7, 2]), variables))

        filter = filters.IndividualDegreeFilter("x-*-A+B", "2-")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'a', 'B2', 'a', 'x-intra-A+B'),
                     Variable('B1', 'a', 'B2', 'a', 'x-intra-B+B'),
                     Variable('A1', 'a', 'B1', 'b', 'x-inter-A+B')]

        self.assertTrue(filter.keep(Monomial([3, 3, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 3, 1, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 3, 3, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 5, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 3, 0, 8]), variables))

        self.assertFalse(filter.keep(Monomial([2, 3, 0, 3]), variables))
        self.assertFalse(filter.keep(Monomial([3, 2, 0, 3]), variables))
        self.assertFalse(filter.keep(Monomial([3, 3, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([3, 3, 3, 2]), variables))
        self.assertFalse(filter.keep(Monomial([2, 2, 3, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 3, 3]), variables))
        self.assertFalse(filter.keep(Monomial([3, 3, 3, 1]), variables))

    def test_sum_degree_filter(self):
        filter = filters.SumDegreeFilter("*", "2+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+B", "1/3+")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'a', 'B2', 'a', 'x-intra-A+B'),
                     Variable('B1', 'a', 'B2', 'a', 'x-intra-B+B')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 5]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 5]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 2, 3]), variables))
        self.assertFalse(filter.keep(Monomial([2, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([5, 2, 0]), variables))

    def test_degree_filter(self):
        filter = filters.DegreeFilter("*", "2+", "3")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([15, 2, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))

        filter = filters.DegreeFilter("x-intra-A+B", "1/3+", "2-/5+")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'a', 'B2', 'a', 'x-intra-A+B'),
                     Variable('B1', 'a', 'B2', 'a', 'x-intra-B+B')]

        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([4, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([2, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([2, 2, 2]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 2, 3]), variables))
        self.assertFalse(filter.keep(Monomial([2, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([3, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 7, 2]), variables))

    def test_num_fragments_filter(self):
        filter = filters.NumFragmentsFilter("*", "1-/3+")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'b', 'A3', 'b', 'x-intra-A+A'),
                     Variable('A2', 'c', 'A3', 'd', 'x-intra-A+A')]

        self.assertTrue(filter.keep(Monomial([3, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 3]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))

        filter = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]


        self.assertTrue(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 2, 3]), variables))
        self.assertTrue(filter.keep(Monomial([2, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 2]), variables))
        self.assertTrue(filter.keep(Monomial([5, 2, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([5, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 7]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))

    def test_not_filter(self):
        filter = filters.IndividualDegreeFilter("*", "2+")
        not_filter = filters.NotFilter(filter)
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        filter = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        not_filter = filters.NotFilter(filter)
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

    def test_and_filter(self):
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        and_filter = filters.AndFilter(filter1, filter2)
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) or filter2.keep(monomial, variables),
                            and_filter.keep(monomial, variables))

    def test_or_filter(self):
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B", "1-")
        or_filter = filters.OrFilter(filter1, filter2)
        variables = [Variable('A1', 'a', 'B1', 'a', 'x-intra-A+B'),
                     Variable('A1', 'b', 'B2', 'b', 'x-intra-A+B'),
                     Variable('B1', 'c', 'B2', 'd', 'x-intra-B+B')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) and filter2.keep(monomial, variables),
                             or_filter.keep(monomial, variables))

    def test_parse_parens(self):
        filter = filters.parse_filter("(", "sum-degree", "x-intra-*+*", "2+", ")")
        variables = [Variable('A1', 'a', 'A2', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A2', 'a', 'A3', 'a', 'x-intra-A+A'),
                     Variable('A1', 'a', 'A1', 'b', 'x-inter-A+A')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0, 5]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1, 2]), variables))

        self.assertFalse(filter.keep(Monomial([1, 1, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0, 1]), variables))


suite = unittest.TestLoader().loadTestsFromTestCase(TestFilters)