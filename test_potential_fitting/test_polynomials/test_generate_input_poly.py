import unittest
import os

from potential_fitting.polynomials import generate_input_poly
from potential_fitting.exceptions import InconsistentValueError

class TestGenerateInputPoly(unittest.TestCase):

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
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)

        self.assertEqual(len(lines), 4)

    def test_A1B2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)

        self.assertEqual(len(lines), 4)

    def test_A1B2_partly_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

    def test_A1B2_purely_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2.in"))

    def test_A1B2X2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)

        self.assertEqual(len(lines), 4)

    def test_A1B2X2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertEqual(len(lines), 4)

    def test_A1B2X2_partly_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

    def test_A1B2X2_purely_inter(self):
        with self.assertRaises(InconsistentValueError):
            generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2X2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2.in"))

    def test_A1B2_C1D2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D1", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D2", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("D1", "b", "D2", "b", "x-intra-D+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "C1", "b", "x-inter-A+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D1", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D2", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D2", "b", "x-inter-B+D"), lines)

        self.assertEqual(len(lines), 17)

    def test_A1B2_C1D2_none(self):
        generate_input_poly(TestGenerateInputPoly.none_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D1", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D2", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("D1", "b", "D2", "b", "x-intra-D+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "C1", "b", "x-inter-A+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D1", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D2", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_filter['degree', 'x-intra-*+*', '1+', '*']", lines)

        self.assertEqual(len(lines), 18)

    def test_A1B2_C1D2_partly_inter(self):
        generate_input_poly(TestGenerateInputPoly.partly_inter_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D1", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D2", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("D1", "b", "D2", "b", "x-intra-D+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "C1", "b", "x-inter-A+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D1", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D2", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_filter['not', 'degree', 'x-inter-*+*', '1+', '*']", lines)

        self.assertEqual(len(lines), 18)

    def test_A1B2_C1D2_purely_inter(self):
        generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D1", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D2", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("D1", "b", "D2", "b", "x-intra-D+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "C1", "b", "x-inter-A+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D1", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D2", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_filter['degree', 'x-intra-*+*', '1+', '*']", lines)

        self.assertEqual(len(lines), 18)

    def test_A1B2X2_C1D2_all(self):
        generate_input_poly(TestGenerateInputPoly.all_settings, "A1B2X2_C1D2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2_C1D2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2X2_C1D2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertIn("add_molecule['{}']\n".format("A1B2X2"), lines)
        self.assertIn("add_molecule['{}']\n".format("C1D2"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D1", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("C1", "b", "D2", "b", "x-intra-C+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("D1", "b", "D2", "b", "x-intra-D+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "C1", "b", "x-inter-A+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D1", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "D2", "b", "x-inter-A+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "C1", "b", "x-inter-B+C"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D1", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "D2", "b", "x-inter-B+D"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X1", "a", "C1", "b", "x-inter-C+X"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X1", "a", "D1", "b", "x-inter-D+X"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X1", "a", "D2", "b", "x-inter-D+X"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X2", "a", "C1", "b", "x-inter-C+X"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X2", "a", "D1", "b", "x-inter-D+X"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("X2", "a", "D2", "b", "x-inter-D+X"), lines)

        self.assertEqual(len(lines), 23)

    def test_A1B2_A1B2_purely_inter(self):
        generate_input_poly(TestGenerateInputPoly.purely_inter_settings, "A1B2_A1B2", os.path.join(TestGenerateInputPoly.output_dir, "A1B2_A1B2.in"))

        with open(os.path.join(TestGenerateInputPoly.output_dir, "A1B2_A1B2.in"), "r") as in_file:
            lines = [line for line in in_file.readlines() if line != "" and line != "\n"]

        self.assertEqual(lines.count("add_molecule['{}']\n".format("A1B2")), 2)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "a", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "a", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "b", "B1", "b", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "b", "B2", "b", "x-intra-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "b", "B2", "b", "x-intra-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "A1", "b", "x-inter-A+A"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B1", "b", "x-inter-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("A1", "a", "B2", "b", "x-inter-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "A1", "b", "x-inter-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B1", "b", "x-inter-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B1", "a", "B2", "b", "x-inter-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "A1", "b", "x-inter-A+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "B1", "b", "x-inter-B+B"), lines)
        self.assertIn("add_variable['{}', '{}', '{}', '{}', '{}']\n".format("B2", "a", "B2", "b", "x-inter-B+B"), lines)
        self.assertIn("add_filter['degree', 'x-intra-*+*', '1+', '*']", lines)

        self.assertEqual(len(lines), 18)
suite = unittest.TestLoader().loadTestsFromTestCase(TestGenerateInputPoly)