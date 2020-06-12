import unittest, random, os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.polynomials import filters, Variable, Monomial
from potential_fitting.exceptions import FilterBadSyntaxError

class TestFilters(TestCaseWithId):

    def test_parse_no_arguments(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Test for when you give parse_filter() no arguments

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter()

        self.test_passed = True

    def test_degree_wrong_num_args(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for when you give 'degree' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', '*')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', '*', '*')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('degree', '*', '*', '*', '*')

        self.test_passed = True

    def test_ind_degree_wrong_num_args(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for when you give 'ind-degre' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', '*')

        self.test_passed = True

    def test_sum_degree_wrong_num_args(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for when you give 'sum-degree' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree', '*')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('sum-degree', '*', '*', '*')

        self.test_passed = True

    def test_num_fragments_wrong_num_args(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for when you give 'num-fragments' the wrong number of arguments.

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments', '*')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('num-fragments', '*', '*', '*')

        self.test_passed = True

    def test_unrecognized_filter_name(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Test for unrecognized filter name

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('some-bad-filter')

        self.test_passed = True

    def test_dangling_conjunction(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for dangling conjunction

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'and')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'or')

        self.test_passed = True

    def test_parse_dangling_not(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Tests for danging not
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not')

        self.test_passed = True

    def test_parse_bad_syntax_in_parens(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Test for when bad syntax appears within parenthesis

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'ind-degree', '*', '*', 'and', ')')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'ind-degree', '*', ')')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('(', 'some-bad-filter', ')')

        self.test_passed = True

    def test_parse_bad_syntax_in_not(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Test for when bad syntax appears within a not

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'ind-degree', '*')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'not')

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('not', 'some-bad-filter')

        self.test_passed = True

    def test_parse_bad_syntax_after_conjunction(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # Test for when bad syntax appears after a conjunction

        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'and', 'sum-degree')
        with self.assertRaises(FilterBadSyntaxError):
            filters.parse_filter('ind-degree', '*', '*', 'and', 'some-bad-filter')

        self.test_passed = True

    def test_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]
        with self.assertRaises(NotImplementedError):
            filters.Filter().keep(Monomial([0, 1, 1]), variables)

        self.test_passed = True

    def test_parse_individual_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.parse_filter('ind-degree', "*", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        self.test_passed = True

    def test_parse_sum_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.parse_filter("sum-degree", "*", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        self.test_passed = True

    def test_parse_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.parse_filter("degree", "*", "2+", "3")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([15, 2, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))

        self.test_passed = True

    def test_parse_num_fragments_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.parse_filter("num-fragments", "*", "1-/3+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'b', 'A', 3, 'b', 'x-intra-A+A-1'),
                     Variable('A', 2, 'c', 'A', 3, 'd', 'x-intra-A+A-0')]

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

        self.test_passed = True

    def test_parse_not_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.IndividualDegreeFilter("*", "2+")
        not_filter = filters.parse_filter("not", "ind-degree", "*", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        filter = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        not_filter = filters.parse_filter("not", "num-fragments", "x-intra-A+B-*", "1-")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        self.test_passed = True

    def test_parse_and_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        and_filter = filters.parse_filter("ind-degree", "*", "2+", "and", "num-fragments", "x-intra-A+B-*", "1-")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) or filter2.keep(monomial, variables),
                            and_filter.keep(monomial, variables))

        self.test_passed = True

    def test_parse_or_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        or_filter = filters.parse_filter("ind-degree", "*", "2+", "or", "num-fragments", "x-intra-A+B-*", "1-")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) and filter2.keep(monomial, variables),
                             or_filter.keep(monomial, variables))

        self.test_passed = True

    def test_individual_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.IndividualDegreeFilter("*", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+B-*", "1/3+")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'a', 'B', 2, 'a', 'x-intra-A+B-1'),
                     Variable('B', 1, 'a', 'B', 2, 'a', 'x-intra-B+B-1')]

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

        filter = filters.IndividualDegreeFilter("x-*-A+B-*", "2-")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'a', 'B', 2, 'a', 'x-intra-A+B-1'),
                     Variable('B', 1, 'a', 'B', 2, 'a', 'x-intra-B+B-1'),
                     Variable('A', 1, 'a', 'B', 1, 'b', 'x-inter-A+B-0')]

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

        self.test_passed = True

    def test_sum_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.SumDegreeFilter("*", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([0, 5, 0]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+B-*", "1/3+")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'a', 'B', 2, 'a', 'x-intra-A+B-1'),
                     Variable('B', 1, 'a', 'B', 2, 'a', 'x-intra-B+B-1')]

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

        self.test_passed = True

    def test_degree_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.DegreeFilter("*", "2+", "3")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([15, 2, 3]), variables))
        self.assertTrue(filter.keep(Monomial([3, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 1]), variables))

        self.assertFalse(filter.keep(Monomial([2, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 1, 2]), variables))

        filter = filters.DegreeFilter("x-intra-A+B-*", "1/3+", "2-/5+")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'a', 'B', 2, 'a', 'x-intra-A+B-1'),
                     Variable('B', 1, 'a', 'B', 2, 'a', 'x-intra-B+B-1')]

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

        self.test_passed = True

    def test_num_fragments_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.NumFragmentsFilter("*", "1-/3+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'b', 'A', 3, 'b', 'x-intra-A+A-1'),
                     Variable('A', 2, 'c', 'A', 3, 'd', 'x-intra-A+A-0')]

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

        filter = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]


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

        self.test_passed = True

    def test_not_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.IndividualDegreeFilter("*", "2+")
        not_filter = filters.NotFilter(filter)
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        filter = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        not_filter = filters.NotFilter(filter)
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(not filter.keep(monomial, variables), not_filter.keep(monomial, variables))

        self.test_passed = True

    def test_and_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        and_filter = filters.AndFilter(filter1, filter2)
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) or filter2.keep(monomial, variables),
                            and_filter.keep(monomial, variables))

        self.test_passed = True

    def test_or_filter(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter1 = filters.IndividualDegreeFilter("*", "2+")
        filter2 = filters.NumFragmentsFilter("x-intra-A+B-*", "1-")
        or_filter = filters.OrFilter(filter1, filter2)
        variables = [Variable('A', 1, 'a', 'B', 1, 'a', 'x-intra-A+B-1'),
                     Variable('A', 1, 'b', 'B', 2, 'b', 'x-intra-A+B-1'),
                     Variable('B', 1, 'c', 'B', 2, 'd', 'x-intra-B+B-0')]

        for i in range(1000):
            monomial = Monomial([random.randint(0, 5), random.randint(0, 5), random.randint(0, 5)])
            self.assertEqual(filter1.keep(monomial, variables) and filter2.keep(monomial, variables),
                             or_filter.keep(monomial, variables))

        self.test_passed = True

    def test_parse_parens(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.parse_filter("(", "sum-degree", "x-intra-*+*-*", "2+", ")")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 2, 'a', 'A', 3, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 1, 'b', 'x-inter-A+A-0')]

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

        self.test_passed = True

    def test_individual_degree_filter_with_levels(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.IndividualDegreeFilter("x-intra-A+A-1", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+A-1-", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+A-0+", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+A-0-1", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.IndividualDegreeFilter("x-intra-A+A-1", "1+")
        variables = [Variable('A', 1, 'aa', 'A', 2, 'aa', 'x-intra-A+A-2'),
                     Variable('A', 1, 'aa', 'A', 3, 'ab', 'x-intra-A+A-1'),
                     Variable('A', 4, 'ab', 'A', 5, 'ab', 'x-intra-A+A-2')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        self.test_passed = True

    def test_sum_degree_filter_with_levels(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.SumDegreeFilter("x-intra-A+A-1", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 2, 0]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+A-1-", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+A-0+", "3+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+A-0-1", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.SumDegreeFilter("x-intra-A+A-1", "1+")
        variables = [Variable('A', 1, 'aa', 'A', 2, 'aa', 'x-intra-A+A-2'),
                     Variable('A', 1, 'aa', 'A', 3, 'ab', 'x-intra-A+A-1'),
                     Variable('A', 4, 'ab', 'A', 5, 'ab', 'x-intra-A+A-2')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        self.test_passed = True

    def test_num_fragments_filter_with_levels(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.NumFragmentsFilter("x-intra-A+A-1", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 2]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 2, 0]), variables))

        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertFalse(filter.keep(Monomial([1, 3, 2]), variables))
        self.assertFalse(filter.keep(Monomial([3, 1, 2]), variables))

        filter = filters.NumFragmentsFilter("x-intra-A+A-1-", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 3]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.NumFragmentsFilter("x-intra-A+A-0+", "3+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.NumFragmentsFilter("x-intra-A+A-0-1", "1+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.NumFragmentsFilter("x-intra-A+A-1-2", "2+")
        variables = [Variable('A', 1, 'aa', 'A', 2, 'aa', 'x-intra-A+A-2'),
                     Variable('A', 1, 'aa', 'A', 3, 'ab', 'x-intra-A+A-1'),
                     Variable('A', 3, 'ab', 'A', 4, 'ab', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 2]), variables))

        self.test_passed = True

    def test_degree_filter_with_levels(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        filter = filters.DegreeFilter("x-intra-A+A-1", "2+", "2")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([3, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 3]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 2]), variables))

        self.assertFalse(filter.keep(Monomial([2, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 2]), variables))


        filter = filters.DegreeFilter("x-intra-A+A-1-", "2+", "2+")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        filter = filters.DegreeFilter("x-intra-A+A-0+", "3+", "4-5")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 1]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 5, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 6]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 2]), variables))
        self.assertTrue(filter.keep(Monomial([2, 1, 2]), variables))
        self.assertTrue(filter.keep(Monomial([3, 3, 0]), variables))

        self.assertFalse(filter.keep(Monomial([5, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 4, 0]), variables))
        self.assertFalse(filter.keep(Monomial([5, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([2, 3, 0]), variables))

        filter = filters.DegreeFilter("x-intra-A+A-0-1", "1+", "2-")
        variables = [Variable('A', 1, 'a', 'A', 2, 'a', 'x-intra-A+A-1'),
                     Variable('A', 1, 'a', 'A', 3, 'b', 'x-intra-A+A-0'),
                     Variable('A', 3, 'b', 'A', 4, 'b', 'x-intra-A+A-1')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 1, 2]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 0, 1]), variables))

        filter = filters.DegreeFilter("x-intra-A+A-1", "1+", "*")
        variables = [Variable('A', 1, 'aa', 'A', 2, 'aa', 'x-intra-A+A-2'),
                     Variable('A', 1, 'aa', 'A', 3, 'ab', 'x-intra-A+A-1'),
                     Variable('A', 4, 'ab', 'A', 5, 'ab', 'x-intra-A+A-2')]

        self.assertTrue(filter.keep(Monomial([0, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 0]), variables))
        self.assertTrue(filter.keep(Monomial([0, 0, 1]), variables))
        self.assertTrue(filter.keep(Monomial([1, 0, 1]), variables))

        self.assertFalse(filter.keep(Monomial([0, 1, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 3, 0]), variables))
        self.assertFalse(filter.keep(Monomial([0, 2, 1]), variables))
        self.assertFalse(filter.keep(Monomial([1, 1, 2]), variables))

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestFilters)
