import sys
sys.path.insert(0, "../../../src/")

import unittest
import os

class TestMonomers(unittest.TestCase):
    def setUp(self):
        self.tests = ["CO2", "N2O4"]
        os.chdir("monomers")

        for test in self.tests:
            directory = test
            settings = test + "/settings.ini"
            database = test + ".db"
            training_set = test + ".xyz"

            if os.system("python ../../../src/database_initializer.py " + settings+ " " + database + " " + directory) >> 8 != 0:
                self.fail()

            if os.system("python ../../../src/database_filler.py " + settings + " " + database + " " + directory) >> 8 != 0:
                self.fail()

            if os.system("python ../../../src/database_reader.py " + settings + " " + database + " " + training_set) >> 8 != 0:
                self.fail()

    def testcompareDataBases(self):
        for test in self.tests:
            if os.system("ndiff -quiet -abserr 1.e-8 expected/" + test + ".xyz " + test + ".xyz") >> 8 != 0:
                self.fail()

    def tearDown(self):
        for test in self.tests:
            os.remove(test + ".db");
            os.remove(test + ".xyz");
        
        os.remove("timer.dat")
        os.system("find . -type d -empty -delete")
        os.chdir("..")

suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(TestMonomers)])
