import os
import configparser
from configparser import NoSectionError, NoOptionError
from potential_fitting.exceptions import ConfigMissingFileError, ConfigMissingSectionError, ConfigMissingPropertyError

class SettingsReader(object):
    """
    Wrapper class for ConfigParser with added functionality
    """

    def __init__(self, file = None):
        """
        Creates a new SettingsReader

        Args:
            file - the .ini file to read from

        Returns:
            A new SettingsReader object
        """

        # create a new ConfigParser
        self.configparser = configparser.ConfigParser(allow_no_value=False)

        if file is not None
            # confirm that the input file exists
            if not os.path.isfile(file):
                raise ConfigMissingFileError(file)

            self.configparser.read(file)

        # save the filepath so that we can include it in error messages
        if file is None:
            self.file = "None"
        else:
            self.file = file
    
    def set(self, section, prop, value):
        self.configparser.set(section, prop, value)

    def write(self, file):
        self.configparser.write(file)

    def get(self, section, prop, default = None):
        """
        Gets the value of a field as a string.

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for

        Return:
            the value of the given property in the given section as a string
        """

        try:
            return self.configparser.get(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file, section, prop) from None
            else:
                return default

    def getboolean(self, section, prop, default = None):
        """
        Gets the value of a field as a boolean.

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for

        Return:
            the value of the given property in the given section as a boolean
        """

        try:
            return self.configparser.getboolean(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file, section, prop) from None
            else:
                return default

    def getint(self, section, prop, default = None):
        """
        Gets the value of a field as an int.

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for

        Return:
            the value of the given property in the given section as an int
        """

        try:
            return self.configparser.getint(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file, section, prop) from None
            else:
                return default

    def getfloat(self, section, prop, default = None):
        """
        Gets the value of a field as a float.

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for

        Return:
            the value of the given property in the given section as a float
        """

        try:
            return self.configparser.getfloat(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file, section, prop) from None
            else:
                return default

    def getlist(self, section, prop, type = str):
        """
        Gets the value of a field as a list of elements of the given type. This list can be multi-dimensional

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for
            type    - the type to cast each element of the list to, default is str

        Return:
            the value of the given property in the given section a list of the given type
        """
        string = self.get(section, prop)

        return parse_array(string, type)

def parse_array(string, type):
    print(string)
    num_open_brackets = 0
    elements = []
    element = ""
    for character in string:
        if num_open_brackets == 1 and character == ",":
            elements.append(element)
            element = ""
        elif character == "[":
            if num_open_brackets != 0:
                element += "["
            num_open_brackets += 1
        elif character == "]":
            num_open_brackets -= 1
            if num_open_brackets != 0:
                element += "]"
        elif not character == " ":
            element += character

    if num_open_brackets == 0:
        if element is not "":
            elements.append(element)
    else:
        print("Something went wrong while parsing a string into a list")
    print(elements)
    return [parse_array(element, type) if "," in element or "[" in element or "]" in element else type(element) for element in elements]
