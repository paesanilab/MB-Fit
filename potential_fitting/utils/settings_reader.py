# external package imports
import os, configparser
from configparser import NoSectionError, NoOptionError

# absolute module imports
from potential_fitting.exceptions import ConfigMissingFileError, ConfigMissingSectionError, ConfigMissingPropertyError

# local module imports
from . import files

class SettingsReader(object):
    """
    Wrapper class for ConfigParser with added functionality.
    """

    def __init__(self, file_path = None):
        """
        Creates a new SettingsReader.

        Args:
            file_path       - Local path to the ".ini" file to read into this SettingsReader.

        Returns:
            A new SettingsReader object.
        """

        # create a new ConfigParser
        self.configparser = configparser.ConfigParser(allow_no_value=False)

        if file_path is not None:
            # confirm that the input file exists
            if not os.path.isfile(file_path):
                raise ConfigMissingFileError(file_path)

            self.configparser.read(file_path)

        # save the filepath so that we can include it in error messages
        if file_path is None:
            self.file_path = "None"
        else:
            self.file_path = file_path
    
    def set(self, section, prop, value):
        """
        Sets the value of a property in a section.

        Args:
            section         - The section of the property to set.
            prop            - The property to set.
            value           - The value to set the property to.

        Returns:
            None.
        """

        self.configparser.set(section, prop, value)

    def write(self, file_path):
        """
        Writes the sections and properties of this SettingsReader to an output file in the ".ini" format.

        Args:
            file_path       - Local path to the ".ini" file to write the settings file to.

        Returns:
            None.
        """

        self.configparser.write(files.init_file(file))

    def get(self, section, prop, default = None):
        """
        Gets the value of a field as a string.
        
        If default is None, then this will throw a ConfigMissingSectionError if the section does not exist and a
        ConfigMissingPropertyError if the property does not exist in the section.

        Args:
            section         - The section of the settings file to look in.
            prop            - The property of the section to look for.

        Return:
            The value of the given property in the given section as a string, or default if the section or property is
            undefined and default is not None.
        """

        try:
            return self.configparser.get(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file_path, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file_path, section, prop) from None
            else:
                return default

    def getboolean(self, section, prop, default = None):
        """
        Gets the value of a field as a boolean.

        If default is None, then this will throw a ConfigMissingSectionError if the section does not exist and a
        ConfigMissingPropertyError if the property does not exist in the section.

        Args:
            section         - The section of the settings file to look in.
            prop            - The property of the section to look for.

        Return:
            The value of the given property in the given section as a boolean, or default if the section or property is
            undefined and default is not None.
        """

        try:
            return self.configparser.getboolean(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file_path, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file_path, section, prop) from None
            else:
                return default

    def getint(self, section, prop, default = None):
        """
        Gets the value of a field as an int.

        If default is None, then this will throw a ConfigMissingSectionError if the section does not exist and a
        ConfigMissingPropertyError if the prop does not exist in the section.

        Args:
            section         - The section of the settings file to look in.
            prop            - The property of the section to look for.

        Return:
            The value of the given property in the given section as an int, or default if the section or property is
            undefined and default is not None.
        """

        try:
            return self.configparser.getint(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file_path, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file_path, section, prop) from None
            else:
                return default

    def getfloat(self, section, prop, default = None):
        """
        Gets the value of a field as a float.

        If default is None, then this will throw a ConfigMissingSectionError if the section does not exist and a
        ConfigMissingPropertyError if the prop does not exist in the section.

        Args:
            section         - The section of the settings file to look in.
            prop            - The property of the section to look for.

        Return:
            The value of the given property in the given section as a float, or default if the section or property is
            undefined and default is not None.
        """

        try:
            return self.configparser.getfloat(section, prop)
        except NoSectionError:
            if default is None:
                raise ConfigMissingSectionError(self.file_path, section, prop) from None
            else:
                return default
        except NoOptionError:
            if default is None:
                raise ConfigMissingPropertyError(self.file_path, section, prop) from None
            else:
                return default

    def getlist(self, section, prop, type = str, default = None):
        """
        Gets the value of a field as a list of elements of the given type. This list can be multi-dimensional.

        If default is None, then this will throw a ConfigMissingSectionError if the section does not exist and a
        ConfigMissingPropertyError if the prop does not exist in the section.

        Args:
            section         - The section of the settings file to look in.
            prop            - The property of the section to look for.
            type            - The type to cast each element of the list to, default is str.

        Return:
            The value of the given property in the given section as a possibly multi-dimensional list of the given
            type, or default if the section or property is undefined and default is not None.
        """

        try:
            string = self.get(section, prop)
        except ConfigMissingSectionError:
            if default is None:
                raise
            else:
                return default
        except NoOptionError:
            if default is None:
                raise
            else:
                return default

        return parse_array(string, type)

def parse_array(string, type):
    """
    Transforms a string into a multidimensional list of the given type.

    Args:
        string              - The string to parse into a list.
        type                - The type of the elements in the list.

    Returns:
        A possibly multi-dimensional list of the given type.
    """

    # keeps track of the number of open braces encountered minus the number of closed braces
    num_open_brackets = 0
    # will be filled with elements
    elements = []
    # characters of the current element
    element = ""

    for character in string:

        # if there is only 1 open brace, and this character is a comma, we reached the end of this element
        if num_open_brackets == 1 and character == ",":
            elements.append(element)
            element = ""

        # if the character is an open brace, and is not the first one, then add it to the element
        elif character == "[":
            if num_open_brackets != 0:
                element += "["
            num_open_brackets += 1

        # if the character is a close brace and it is not the last one, then add it to the element
        elif character == "]":
            num_open_brackets -= 1
            if num_open_brackets != 0:
                element += "]"

        # if the character is not whitespace, add it to the current element
        elif not character == " ":
            element += character

    # if we closed all opened braces and the last element isn't empty, add it to the list of elements
    if num_open_brackets == 0:
        if element is not "":
            elements.append(element)

    # if we didn't close all braces, there is a problem
    else:
        print("Something went wrong while parsing a string into a list")

    # perfrom recursive call on all elements that still have a ",", "[" or "]", otherwise, just cast to the given type
    return [parse_array(element, type) if "," in element or "[" in element or "]" in element else type(element) for element in elements]
