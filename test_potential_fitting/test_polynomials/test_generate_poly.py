import unittest
import os

from potential_fitting.polynomials import generate_poly

class TestGeneratePoly(unittest.TestCase):

    def setUpClass():
        TestGeneratePoly.resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
        TestGeneratePoly.output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
        TestGeneratePoly.reference_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
        TestGeneratePoly.settings_path = os.path.join(TestGeneratePoly.resources_dir, "none.ini")

    def test_A4B1_degree4_no_filter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A4B1_no_filter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A4B1_no_filter_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A4B1_no_filter_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1B3_degree4_no_filter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B3_no_filter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B3_no_filter_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B3_no_filter_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A2B1_A2B1_degree2_partly_inter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A2B1_A2B1_partly_inter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A2B1_A2B1_partly_inter_degree_2")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A2B1_A2B1_partly_inter_degree_2")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      2,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A2B5_C1D2_degree2_partly_inter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A2B5_C1D2_partly_inter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A2B5_C1D2_partly_inter_degree_2")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A2B5_C1D2_partly_inter_degree_2")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      2,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1B2X2_A1B2X2_degree4_no_filters(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B2X2_A1B2X2_no_filters.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B2X2_A1B2X2_no_filters_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B2X2_A1B2X2_no_filters_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1_B1C2_B1C2_degree2_partly_inter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1_B1C2_B1C2_partly_inter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1_B1C2_B1C2_partly_inter_degree_2")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1_B1C2_B1C2_partly_inter_degree_2")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      2,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())


    def test_A1B2X2_A1B2X2_A1B2X2_degree2_purely_inter(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B2X2_A1B2X2_A1B2X2_purely_inter.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B2X2_A1B2X2_A1B2X2_purely_inter_degree_2")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B2X2_A1B2X2_A1B2X2_purely_inter_degree_2")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      2,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1B2C2_degree4_custom_filters(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B2C2_custom_filters.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B2C2_custom_filters_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B2C2_custom_filters_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1B2X2_A1B2X2_degree4_mbpol_filters(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B2X2_A1B2X2_mbpol_filters.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B2X2_A1B2X2_mbpol_filters_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B2X2_A1B2X2_mbpol_filters_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())

    def test_A1B2_A1B2_A1B2_degree4_mbpol_filters(self):

        input_path = os.path.join(TestGeneratePoly.resources_dir, "A1B2_A1B2_A1B2_mbpol_filters.in")
        output_dir = os.path.join(TestGeneratePoly.output_dir, "A1B2_A1B2_A1B2_mbpol_filters_degree_4")
        reference_dir = os.path.join(TestGeneratePoly.reference_dir, "A1B2_A1B2_A1B2_mbpol_filters_degree_4")

        generate_poly(TestGeneratePoly.settings_path,
                      input_path,
                      4,
                      output_dir
                      )

        out_log_path = os.path.join(output_dir, "poly.log")
        ref_log_path = os.path.join(reference_dir, "poly.log")

        with open(out_log_path, "r") as out_log, open(ref_log_path, "r") as ref_log:
            self.assertEqual(out_log.readlines(), ref_log.readlines())



suite = unittest.TestLoader().loadTestsFromTestCase(TestGeneratePoly)