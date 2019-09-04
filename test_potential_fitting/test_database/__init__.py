import unittest
from . import test_database, test_database_job_reader_and_writer

suite = unittest.TestSuite([test_database.suite, test_database_job_reader_and_writer.suite])
