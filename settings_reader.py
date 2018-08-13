import configparser
from configparser import NoSectionError, NoOptionError
from exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

class SettingsReader(object):
    
    def __init__(self, file):
        self.configparser = configparser.ConfigParser(allow_no_value=False)
        self.configparser.read(file)

        self.file = file

    def get(self, section, prop):
        try:
            return self.configparser.get(section, prop)
        except NoSectionError:
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None

    def getboolean(self, section, prop):
        try:
            return self.configparser.getboolean(section, prop)
        except NoSectionError:
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None

    def getint(self, section, prop):
        try:
            return self.configparser.getint(section, prop)
        except NoSectionError:
            raise ConfigMissingSectionError(self.file, section, prop) from None
        except NoOptionError:
            raise ConfigMissingPropertyError(self.file, section, prop) from None
