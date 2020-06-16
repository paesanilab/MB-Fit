import unittest, os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.utils import files
from potential_fitting.exceptions import FileExistsError, InvalidValueError

class TestFiles(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestFiles, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    def setUpClass():
        TestFiles.dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dir")
        TestFiles.dirdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dir", "dir")

        TestFiles.dir2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dir2")
        TestFiles.dir2dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "dir2", "dir")

        TestFiles.file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "file1.txt")
        TestFiles.file_backup1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "file1.txt.backup-1")
        TestFiles.file_backup2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "file1.txt.backup-2")
        TestFiles.file_backup3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "file1.txt.backup-3")


        TestFiles.deepdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "some", "file", "path")

        TestFiles.deepfile = os.path.join(TestFiles.deepdir, "file.txt")

    def test_init_directory(self):

        self.assertEqual(files.init_directory(TestFiles.dir), TestFiles.dir)
        self.assertTrue(os.path.isdir(TestFiles.dir))
        self.assertEqual(files.init_directory(TestFiles.dir), TestFiles.dir)
        self.assertTrue(os.path.isdir(TestFiles.dir))


        self.assertEqual(files.init_directory(TestFiles.dirdir), TestFiles.dirdir)
        self.assertTrue(os.path.isdir(TestFiles.dirdir))

        self.assertEqual(files.init_directory(TestFiles.dir2dir), TestFiles.dir2dir)
        self.assertTrue(os.path.isdir(TestFiles.dir2dir))

        os.removedirs(TestFiles.dirdir)
        os.removedirs(TestFiles.dir2dir)

        self.test_passed = True

    def test_init_file(self):
        self.assertEqual(files.init_file(TestFiles.file), TestFiles.file)
        self.assertFalse(os.path.isfile(TestFiles.file))

        with open(TestFiles.file, "w") as file:
            file.write("FILE #1.")

        self.assertEqual(files.init_file(TestFiles.file), TestFiles.file)
        self.assertFalse(os.path.isfile(TestFiles.file))
        self.assertTrue(os.path.isfile(TestFiles.file_backup1))

        with open(TestFiles.file_backup1, "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open(TestFiles.file, "w") as file:
            file.write("FILE #2.")

        self.assertEqual(files.init_file(TestFiles.file), TestFiles.file)
        self.assertFalse(os.path.isfile(TestFiles.file))
        self.assertTrue(os.path.isfile(TestFiles.file_backup1))
        self.assertTrue(os.path.isfile(TestFiles.file_backup2))

        with open(TestFiles.file_backup1, "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open(TestFiles.file_backup2, "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        with open(TestFiles.file, "w") as file:
            file.write("FILE #3.")

        self.assertEqual(files.init_file(TestFiles.file, files.OverwriteMethod.OVERWRITE), TestFiles.file)
        self.assertTrue(os.path.isfile(TestFiles.file))
        self.assertTrue(os.path.isfile(TestFiles.file_backup1))
        self.assertTrue(os.path.isfile(TestFiles.file_backup2))

        with open(TestFiles.file, "r") as file:
            self.assertEqual(file.read(), "FILE #3.")

        with open(TestFiles.file_backup1, "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open(TestFiles.file_backup2, "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        with open(TestFiles.file, "w") as file:
            file.write("FILE #4.")

        with self.assertRaises(FileExistsError):
            files.init_file(TestFiles.file, files.OverwriteMethod.CRASH)

        with open(TestFiles.file, "r") as file:
            self.assertEqual(file.read(), "FILE #4.")

        with open(TestFiles.file_backup1, "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open(TestFiles.file_backup2, "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        self.assertEqual(files.init_file(TestFiles.file, files.OverwriteMethod.NONE), TestFiles.file)
        self.assertTrue(os.path.isfile(TestFiles.file))
        self.assertTrue(os.path.isfile(TestFiles.file_backup1))
        self.assertTrue(os.path.isfile(TestFiles.file_backup2))

        with open(TestFiles.file, "r") as file:
            self.assertEqual(file.read(), "FILE #4.")

        with open(TestFiles.file_backup1, "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open(TestFiles.file_backup2, "r") as file:
            self.assertEqual(file.read(), "FILE #2.")


        self.assertEqual(files.init_file(TestFiles.deepfile), TestFiles.deepfile)
        self.assertTrue(os.path.isdir(TestFiles.deepdir))

        with self.assertRaises(InvalidValueError):
            files.init_file(TestFiles.file, 8)

        os.remove(TestFiles.file)
        os.remove(TestFiles.file_backup1)
        os.remove(TestFiles.file_backup2)
        os.removedirs(TestFiles.deepdir)

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestFiles)
