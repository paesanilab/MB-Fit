import configparser
from configparser import NoSectionError, NoOptionError
from exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

class SettingsReader(object):
    """
    Wrapper class for ConfigParser with added functionality
    """

    def __init__(self, file):
        """
        Creates a new SettingsReader

        Args:
            file - the .ini file to read from

        Returns:
            A new SettingsReader object
        """

        self.configparser = configparser.ConfigParser(allow_no_value=False)
        self.configparser.read(file)

        # save the filepath so that we can include it in error messages
        self.file = file

    def get(self, section, prop):
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
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None

    def getboolean(self, section, prop):
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
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None

    def getint(self, section, prop):
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
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None

    def getlist(self, section, prop):
        """
        Gets the value of a field as a list of strings. This list can be multi-dimensional

        Will throw a ConfigMissingSectionError if the section does not exist
        and a ConfigMissingPropertyError if the prop does not exist in the section

        Args:
            section - the section of the settings file to look in
            prop    - the property of the section to look for

        Return:
            the value of the given property in the given section a list of strings.
        """
        string = self.get(section, prop)
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
            else:
                element += character

        if num_open_brackets == 0:
            elements.append(element)
        else:
            print("Something went wrong while parsing a string into a list")
        return [parse_array(element) if "," in element else element for element in elements]
