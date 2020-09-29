import unittest
import os

from test_mbfit.test_case_with_id import TestCaseWithId
from mbfit.polynomials import generate_input_poly
from mbfit.exceptions import InconsistentValueError

class TestGenerateInputPoly(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestGenerateInputPoly, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    def setUpClass():
        TestGenerateInputPoly.all_settings = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "all.ini")
        TestGenerateInputPoly.none_settings = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "none.ini")
        TestGenerateInputPoly.partly_inter_settings = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "partly-inter.ini")
        TestGenerateInputPoly.purely_inter_settings = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "purely-inter.ini")

        TestGenerateInputPoly.output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

    def test_A1B2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)

        self.assertEqual(len(lines), 4)

        self.test_passed = True

    def test_A1B2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)

        self.assertEqual(len(lines), 4)

        self.test_passed = True

    def test_A1B2_partly_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

        self.test_passed = True

    def test_A1B2_purely_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

        self.test_passed = True

    def test_A1B2X2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)

        self.assertEqual(len(lines), 4)

        self.test_passed = True

    def test_A1B2X2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertEqual(len(lines), 4)

        self.test_passed = True

    def test_A1B2X2_partly_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        self.test_passed = True

    def test_A1B2X2_purely_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        self.test_passed = True

    def test_A1B2_C1D2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 1, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 2, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "D", 2, "b", "x-intra-D+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "C", 1, "b", "x-inter-A+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 1, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 2, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 2, "b", "x-inter-B+D-0"), lines)

        self.assertEqual(len(lines), 17)

        self.test_passed = True

    def test_A1B2_C1D2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 1, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 2, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "D", 2, "b", "x-intra-D+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "C", 1, "b", "x-inter-A+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 1, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 2, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_filter['sum-degree', 'x-inter-*+*-*', '0']", lines)
        self.assertEqual(len(lines), 18)

        self.test_passed = True

    def test_A1B2_C1D2_partly_inter(self):
        generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 1, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 2, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "D", 2, "b", "x-intra-D+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "C", 1, "b", "x-inter-A+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 1, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 2, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_filter['sum-degree', 'x-inter-*+*-*', '0']", lines)

        self.assertEqual(len(lines), 18)

        self.test_passed = True

    def test_A1B2_C1D2_purely_inter(self):
        generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 1, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 2, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "D", 2, "b", "x-intra-D+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "C", 1, "b", "x-inter-A+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 1, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 2, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_filter['ind-degree', 'x-intra-*+*-*', '1+']", lines)

        self.assertEqual(len(lines), 18)

        self.test_passed = True

    def test_A1B2X2_C1D2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2X2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 1, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "D", 2, "b", "x-intra-C+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "D", 2, "b", "x-intra-D+D-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "C", 1, "b", "x-inter-A+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 1, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "D", 2, "b", "x-inter-A+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "C", 1, "b", "x-inter-B+C-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 1, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "D", 2, "b", "x-inter-B+D-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "X", 1, "a", "x-inter-C+X-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("C", 1, "b", "X", 2, "a", "x-inter-C+X-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "X", 1, "a", "x-inter-D+X-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 1, "b", "X", 2, "a", "x-inter-D+X-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 2, "b", "X", 1, "a", "x-inter-D+X-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("D", 2, "b", "X", 2, "a", "x-inter-D+X-0"), lines)

        self.assertEqual(len(lines), 23)

        self.test_passed = True

    def test_A1B2_A1B2_purely_inter(self):
        generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2_A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_A1B2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_A1B2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertEqual(lines.count("add_molecule['{}']\n".format("A1B2")), 2)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 1, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 2, "a", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 2, "a", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 2, "b", "B", 3, "b", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 2, "b", "B", 4, "b", "x-intra-A+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 3, "b", "B", 4, "b", "x-intra-B+B-1"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "A", 2, "b", "x-inter-A+A-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 3, "b", "x-inter-A+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 1, "a", "B", 4, "b", "x-inter-A+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 2, "b", "B", 1, "a", "x-inter-A+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 3, "b", "x-inter-B+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 1, "a", "B", 4, "b", "x-inter-B+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("A", 2, "b", "B", 2, "a", "x-inter-A+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "B", 3, "b", "x-inter-B+B-0"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format("B", 2, "a", "B", 4, "b", "x-inter-B+B-0"), lines)
        self.assertIn("add_filter['ind-degree', 'x-intra-*+*-*', '1+']", lines)

        self.assertEqual(len(lines), 18)

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestGenerateInputPoly)
