import unittest, os

from test_mbfit.test_case_with_id import TestCaseWithId
from mbfit.utils import settings_reader
from mbfit.exceptions import ConfigMissingFileError, ConfigMissingSectionError, ConfigMissingPropertyError

class TestSettingsReader(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestSettingsReader, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        self.settings_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "settings.ini")
        self.settings_reader = settings_reader.SettingsReader(self.settings_file)

    def test_constructor_and_get_file_path(self):

        with self.assertRaises(ConfigMissingFileError):
            settings_reader.SettingsReader(
                os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "no_file.ini"))

        self.assertIsNone(settings_reader.SettingsReader().get_file_path())
        self.assertEqual(self.settings_reader.get_file_path(), self.settings_file)

        self.test_passed = True

    def test_get(self):

        self.assertEqual(self.settings_reader.get("section1", "property1"), "value1")

        with self.assertRaises(ConfigMissingSectionError):
            self.settings_reader.get("no_section", "no_property")

        self.assertEqual(self.settings_reader.get("no_section", "no_property", "default"), "default")

        with self.assertRaises(ConfigMissingPropertyError):
            self.settings_reader.get("section1", "no_property")

        self.assertEqual(self.settings_reader.get("section1", "no_property", "default"), "default")

        self.test_passed = True

    def test_getboolean(self):

        self.assertEqual(self.settings_reader.getboolean("section1", "boolean1"), True)
        self.assertEqual(self.settings_reader.getboolean("section1", "boolean2"), False)

        with self.assertRaises(ConfigMissingSectionError):
            self.settings_reader.getboolean("no_section", "no_property")

        self.assertEqual(self.settings_reader.getboolean("no_section", "no_property", True), True)

        with self.assertRaises(ConfigMissingPropertyError):
            self.settings_reader.getboolean("section1", "no_property")

        self.assertEqual(self.settings_reader.getboolean("section1", "no_property", False), False)

        self.test_passed = True

    def test_getint(self):

        self.assertEqual(self.settings_reader.getint("section1", "int"), 28)

        with self.assertRaises(ConfigMissingSectionError):
            self.settings_reader.getint("no_section", "no_property")

        self.assertEqual(self.settings_reader.getint("no_section", "no_property", 14), 14)

        with self.assertRaises(ConfigMissingPropertyError):
            self.settings_reader.getint("section1", "no_property")

        self.assertEqual(self.settings_reader.getint("section1", "no_property", 98), 98)

        self.test_passed = True

    def test_getfloat(self):

        self.assertEqual(self.settings_reader.getfloat("section1", "float"), 6.626)

        with self.assertRaises(ConfigMissingSectionError):
            self.settings_reader.getfloat("no_section", "no_property")

        self.assertEqual(self.settings_reader.getfloat("no_section", "no_property", 14.32), 14.32)

        with self.assertRaises(ConfigMissingPropertyError):
            self.settings_reader.getfloat("section1", "no_property")

        self.assertEqual(self.settings_reader.getfloat("section1", "no_property", 98.96), 98.96)

        self.test_passed = True

    def test_getlist(self):
        self.assertEqual(self.settings_reader.getlist("section1", "str_list"), [["red", "fish"], ["blue", "fish"]])
        self.assertEqual(self.settings_reader.getlist("section1", "int_list", int), [1, 4, 3, 2])

        with self.assertRaises(ConfigMissingSectionError):
            self.settings_reader.getlist("no_section", "no_property")

        self.assertEqual(self.settings_reader.getlist("no_section", "no_property", float, [14.32]), [14.32])

        with self.assertRaises(ConfigMissingPropertyError):
            self.settings_reader.getlist("section1", "no_property")

        self.assertEqual(self.settings_reader.getlist("section1", "no_property", str, ["list", "item"]), ["list", "item"])

        self.test_passed = True

    def test_set(self):
        self.settings_reader.set("section1", "property1", "different_value!")
        self.assertEqual(self.settings_reader.get("section1", "property1"), "different_value!")

        self.settings_reader.set("section1", "property2", "value2")
        self.assertEqual(self.settings_reader.get("section1", "property2"), "value2")

        self.settings_reader.set("section2", "property1", "value3")
        self.assertEqual(self.settings_reader.get("section2", "property1"), "value3")

        self.test_passed = True

    def test_write(self):
        s_r = settings_reader.SettingsReader()
        s_r.set("section1", "property1", "value1")
        s_r.set("section2", "property2", "value2")

        write_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "write.ini")

        s_r.write(write_path)

        s_r = settings_reader.SettingsReader(write_path)

        self.assertEqual(s_r.get("section1", "property1"), "value1")
        self.assertEqual(s_r.get("section2", "property2"), "value2")

        os.remove(os.path.join(os.path.dirname(os.path.abspath(__file__)), "write.ini"))

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestSettingsReader)
