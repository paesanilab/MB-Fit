import unittest, os

from potential_fitting.database import Database


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


suite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)