import unittest, os, psycopg2

from potential_fitting.database import Database
from potential_fitting.exceptions import InvalidValueError


class TestDatabase(unittest.TestCase):

    def setUp(self):
        self.config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database.ini)")
        self.database = Database(self.config)

    def tearDown(self):
        self.database.annihilate(confirm="confirm")
        self.database.close()

    def test_set_and_get_batch_size(self):
        self.database.set_batch_size(8)
        self.assertEqual(self.database.get_batch_size(), 8)

        self.database.set_batch_size(15)
        self.assertEqual(self.database.get_batch_size(), 15)

        with self.assertRaises(InvalidValueError):
            self.database.set_batch_size(-1)

        self.assertEqual(self.database.get_batch_size(), 15)

    def test_create(self):
        pass

    def test_execute(self):

        with self.assertRaises(psycopg2.InternalError):
            self.database.execute("RAISE EXCEPTION 'EXECUTE WORKED'")
        self.database.save()

        pass

    def test_create_postgres_array(self):
        self.assertEqual(self.database.create_postgres_array("A", "B", "C", "D"), "{A, B, C, D}")
        self.assertEqual(self.database.create_postgres_array(1, 2, 3, 4), "{1, 2, 3, 4}")


suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabase)