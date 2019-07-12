import unittest, os

from potential_fitting.utils import files
from potential_fitting.exceptions import FileExistsError

class TestFiles(unittest.TestCase):

    def test_init_directory(self):

        self.assertEqual(files.init_directory("dir"), "dir")
        self.assertTrue(os.path.isdir("dir"))
        self.assertEqual(files.init_directory("dir"), "dir")
        self.assertTrue(os.path.isdir("dir"))


        self.assertEqual(files.init_directory("dir/dir"), "dir/dir")
        self.assertTrue(os.path.isdir("dir/dir"))

        self.assertEqual(files.init_directory("dir2/dir"), "dir2/dir")
        self.assertTrue(os.path.isdir("dir2/dir"))

        os.removedirs("dir/dir")
        os.removedirs("dir2/dir")

    def test_init_file(self):
        self.assertEqual(files.init_file("file.txt"), "./file.txt")
        self.assertFalse(os.path.isfile("file.txt"))

        with open("file.txt", "w") as file:
            file.write("FILE #1.")

        self.assertEqual(files.init_file("./file.txt"), "./file.txt")
        self.assertFalse(os.path.isfile("file.txt"))
        self.assertTrue(os.path.isfile("file.txt.backup-1"))

        with open("file.txt.backup-1", "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open("file.txt", "w") as file:
            file.write("FILE #2.")

        self.assertEqual(files.init_file("./file.txt"), "./file.txt")
        self.assertFalse(os.path.isfile("file.txt"))
        self.assertTrue(os.path.isfile("file.txt.backup-1"))
        self.assertTrue(os.path.isfile("file.txt.backup-2"))

        with open("file.txt.backup-1", "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open("file.txt.backup-2", "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        with open("file.txt", "w") as file:
            file.write("FILE #3.")

        self.assertEqual(files.init_file("file.txt", files.OverwriteMethod.OVERWRITE), "./file.txt")
        self.assertTrue(os.path.isfile("file.txt"))
        self.assertTrue(os.path.isfile("file.txt.backup-1"))
        self.assertTrue(os.path.isfile("file.txt.backup-2"))

        with open("file.txt", "r") as file:
            self.assertEqual(file.read(), "FILE #3.")

        with open("file.txt.backup-1", "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open("file.txt.backup-2", "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        with open("file.txt", "w") as file:
            file.write("FILE #4.")

        with self.assertRaises(FileExistsError):
            files.init_file("file.txt", files.OverwriteMethod.CRASH)

        with open("file.txt", "r") as file:
            self.assertEqual(file.read(), "FILE #4.")

        with open("file.txt.backup-1", "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open("file.txt.backup-2", "r") as file:
            self.assertEqual(file.read(), "FILE #2.")

        self.assertEqual(files.init_file("file.txt", files.OverwriteMethod.NONE), "./file.txt")
        self.assertTrue(os.path.isfile("file.txt"))
        self.assertTrue(os.path.isfile("file.txt.backup-1"))
        self.assertTrue(os.path.isfile("file.txt.backup-2"))

        with open("file.txt", "r") as file:
            self.assertEqual(file.read(), "FILE #4.")

        with open("file.txt.backup-1", "r") as file:
            self.assertEqual(file.read(), "FILE #1.")

        with open("file.txt.backup-2", "r") as file:
            self.assertEqual(file.read(), "FILE #2.")


        self.assertEqual(files.init_file("some/path/file.txt"), "some/path/file.txt")
        self.assertTrue(os.path.isdir("some/path/"))

        os.remove("file.txt")
        os.remove("file.txt.backup-1")
        os.remove("file.txt.backup-2")
        os.removedirs("some/path/")


suite = unittest.TestLoader().loadTestsFromTestCase(TestFiles)