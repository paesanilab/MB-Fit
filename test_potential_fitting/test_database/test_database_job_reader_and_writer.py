
import unittest, os
from glob import glob

from potential_fitting import database
from potential_fitting.exceptions import DatabaseConnectionError
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

    def setUp(self):
        self.config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user1.ini")
        self.water_opt_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_opt.xyz")
        self.dimer_configs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_dimer_configs.xyz")
        self.mon_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water.ini")
        self.dim_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_dimer.ini")
        self.execute_jobs_script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "execute_jobs.sh")
        self.job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "jobs")
        self.dim_train_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dim_train_set.xyz")

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

    def test_job_writer_and_reader(self):
        database.initialize_database(self.mon_settings_path, self.config_path, self.water_opt_path, "HF", "STO-3G", False, "database_test", optimized=True)
        database.initialize_database(self.dim_settings_path, self.config_path, self.dimer_configs_path, "HF", "STO-3G", False, "database_test")

        # running this should make all 10 jobs
        database.make_all_jobs(self.mon_settings_path, self.config_path, "testclient", self.job_dir, "database_test", num_jobs=10)
        self.assertEqual(len(glob(self.job_dir + "/job_*.py")), 10)

        # no more jobs should be made because all 10 have been made
        database.make_all_jobs(self.mon_settings_path, self.config_path, "testclient", self.job_dir, "database_test", num_jobs=10)
        self.assertEqual(len(glob(self.job_dir + "/job_*.py")), 10)

        system.call(self.execute_jobs_script_path, self.job_dir, out_file=None)

        self.assertEqual(len(glob(self.job_dir + "/job_*")), 20)

        database.read_all_jobs(self.config_path, self.job_dir)


        self.assertEqual(len(glob(self.job_dir + "/job_*")), 20)
        self.assertEqual(len(glob(self.job_dir + "/job_*_done")), 10)

        database.generate_training_set(self.dim_settings_path, self.config_path, self.dim_train_path, "HF", "STO-3G", False, "database_test")

        config_molecules = list(parse_training_set_file(self.dimer_configs_path, settings=SettingsReader(self.dim_settings_path)))
        train_molecules = list(parse_training_set_file(self.dim_train_path, settings=SettingsReader(self.dim_settings_path)))

        for molecule in train_molecules:
            self.assertIn(molecule, config_molecules)
        for molecule in config_molecules:
            self.assertIn(molecule, train_molecules)





suite = unittest.TestLoader().loadTestsFromTestCase(TestJobReaderAndWriter)
