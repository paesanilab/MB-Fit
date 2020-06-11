import unittest, os
from glob import glob

from potential_fitting import database
from potential_fitting.exceptions import DatabaseConnectionError, CommandExecutionError
from potential_fitting.utils import system, SettingsReader
from potential_fitting.molecule import parse_training_set_file

# only import psycopg2 if it is installed.
try:
    import psycopg2
except ModuleNotFoundError:
    pass

def psycopg2_installed():
    try:
        import psycopg2
        return True
    except ModuleNotFoundError:
        return False

def hasQchem():

    try:
        system.call("which", "qchem")
    except CommandExecutionError:
        return False
    return True

def hasPsi4():

    try:
        import psi4
    except ImportError:
        return False
    return True

def local_db_installed():
    try:
        config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user1.ini")
        db = database.Database(config)
        db.close()
        return True
    except DatabaseConnectionError:
        return False

@unittest.skipUnless(psycopg2_installed() and local_db_installed(),"psycopg2 or a local test database is not installed, so database cannot be tested.")
class TestJobReaderAndWriter(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestJobReaderAndWriter, self).__init__(*args, **kwargs)
        self.test_passed = False
        self.test_name = self.id()

    def setUpClass():
        TestJobReaderAndWriter.config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user1.ini")
        TestJobReaderAndWriter.water_opt_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_opt.xyz")
        TestJobReaderAndWriter.dimer_configs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_dimer_configs.xyz")
        TestJobReaderAndWriter.mon_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water.ini")
        TestJobReaderAndWriter.psi4_dim_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_dimer_psi4.ini")
        TestJobReaderAndWriter.qchem_dim_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_dimer_qchem.ini")
        TestJobReaderAndWriter.execute_jobs_script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "execute_jobs.sh")
        TestJobReaderAndWriter.job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "jobs")
        TestJobReaderAndWriter.dim_train_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dim_train_set.xyz")

    def setUp(self):

        system.call("rm", self.job_dir, "-r", "-f")

        db = database.Database(self.config_path)
        db.annihilate(confirm="confirm")
        db.save()
        db.close()

    def tearDown(self):
        db = database.Database(self.config_path)
        db.annihilate(confirm="confirm")
        db.save()
        db.close()

        system.call("rm", self.job_dir, "-r", "-f")

        local_output = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output")
        mbfithome = os.environ.get('MBFIT_HOME')
        if os.path.isdir(local_output):
            if self.test_passed:
                os.system("mkdir -p " + os.path.join(mbfithome, "passed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "passed_tests_outputs", self.test_name))
            else:
                os.system("mkdir -p " + os.path.join(mbfithome, "failed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "failed_tests_outputs", self.test_name))

    @unittest.skipUnless(hasQchem(), "Qchem is not installed and cannot be tested.")
    def test_job_writer_and_reader_qchem(self):
        database.initialize_database(TestJobReaderAndWriter.mon_settings_path,
                                     TestJobReaderAndWriter.config_path,
                                     TestJobReaderAndWriter.water_opt_path,
                                     "HF", "STO-3G", False, "database_test", optimized=True)
        database.initialize_database(TestJobReaderAndWriter.qchem_dim_settings_path,
                                     TestJobReaderAndWriter.config_path,
                                     TestJobReaderAndWriter.dimer_configs_path,
                                     "HF", "STO-3G", False, "database_test")

        # running this should make all 10 jobs
        job_handler = database.get_job_handler(TestJobReaderAndWriter.qchem_dim_settings_path)
        job_handler.make_all_jobs(TestJobReaderAndWriter.config_path,
                                  "testclient", self.job_dir, "database_test", num_jobs=10)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*.py")), 10)

        # no more jobs should be made because all 10 have been made
        job_handler.make_all_jobs(TestJobReaderAndWriter.config_path,
                                  "testclient", TestJobReaderAndWriter.job_dir,
                                  "database_test", num_jobs=10)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*.py")), 10)

        system.call(TestJobReaderAndWriter.execute_jobs_script_path, TestJobReaderAndWriter.job_dir, out_file=None)

        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*")), 20)

        job_handler.read_all_jobs(TestJobReaderAndWriter.config_path, TestJobReaderAndWriter.job_dir)

        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*")), 20)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*_done")), 10)

        database.generate_training_set(TestJobReaderAndWriter.qchem_dim_settings_path,
                                       TestJobReaderAndWriter.config_path,
                                       TestJobReaderAndWriter.dim_train_path,
                                       "HF", "STO-3G", False, "database_test")

        config_molecules = list(parse_training_set_file(TestJobReaderAndWriter.dimer_configs_path,
                                                        settings=SettingsReader(TestJobReaderAndWriter.qchem_dim_settings_path)))
        train_molecules = list(parse_training_set_file(TestJobReaderAndWriter.dim_train_path,
                                                       settings=SettingsReader(TestJobReaderAndWriter.qchem_dim_settings_path)))

        for molecule in train_molecules:
            self.assertIn(molecule, config_molecules)
        for molecule in config_molecules:
            self.assertIn(molecule, train_molecules)

        self.test_passed = True

    @unittest.skipUnless(hasPsi4(), "Psi4 is not installed and cannot be tested.")
    def test_job_writer_and_reader_psi4(self):
        database.initialize_database(TestJobReaderAndWriter.mon_settings_path,
                                     TestJobReaderAndWriter.config_path,
                                     TestJobReaderAndWriter.water_opt_path,
                                     "HF", "STO-3G", False, "database_test", optimized=True)
        database.initialize_database(TestJobReaderAndWriter.psi4_dim_settings_path,
                                     TestJobReaderAndWriter.config_path,
                                     TestJobReaderAndWriter.dimer_configs_path,
                                     "HF", "STO-3G", False, "database_test")

        # running this should make all 10 jobs
        job_handler = database.get_job_handler(TestJobReaderAndWriter.psi4_dim_settings_path)
        job_handler.make_all_jobs(TestJobReaderAndWriter.config_path,
                                  "testclient", self.job_dir, "database_test", num_jobs=10)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*.py")), 10)

        # no more jobs should be made because all 10 have been made
        job_handler.make_all_jobs(TestJobReaderAndWriter.config_path,
                                  "testclient", TestJobReaderAndWriter.job_dir,
                                  "database_test", num_jobs=10)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*.py")), 10)

        system.call(TestJobReaderAndWriter.execute_jobs_script_path, TestJobReaderAndWriter.job_dir, out_file=None)

        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*")), 20)

        job_handler.read_all_jobs(TestJobReaderAndWriter.config_path, TestJobReaderAndWriter.job_dir)

        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*")), 20)
        self.assertEqual(len(glob(TestJobReaderAndWriter.job_dir + "/job_*_done")), 10)

        database.generate_training_set(TestJobReaderAndWriter.psi4_dim_settings_path,
                                       TestJobReaderAndWriter.config_path,
                                       TestJobReaderAndWriter.dim_train_path,
                                       "HF", "STO-3G", False, "database_test")

        config_molecules = list(parse_training_set_file(TestJobReaderAndWriter.dimer_configs_path,
                                                        settings=SettingsReader(TestJobReaderAndWriter.qchem_dim_settings_path)))
        train_molecules = list(parse_training_set_file(TestJobReaderAndWriter.dim_train_path,
                                                       settings=SettingsReader(TestJobReaderAndWriter.psi4_dim_settings_path)))

        for molecule in train_molecules:
            self.assertIn(molecule, config_molecules)
        for molecule in config_molecules:
            self.assertIn(molecule, train_molecules)

        self.test_passed = True





suite = unittest.TestLoader().loadTestsFromTestCase(TestJobReaderAndWriter)
