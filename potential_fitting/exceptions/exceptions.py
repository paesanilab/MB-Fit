"""
--------------------------- General Exceptions
"""

class PotentialFittingError(Exception):
    """Basic exception for all errors raised by our code"""
    
    def __init__(self, message):
        super().__init__("The following error occured in the Potential Fitting Library: {}".format(message))


"""
--------------------------- File Errors
"""

class FileExistsError(PotentialFittingError):
    """Raised when a file must be overwritten but it already exists"""

    def __init__(self, file):
        super().__init__("File '{}' already exists, but must be overwritten".format(file))

"""
--------------------------- Called Process Errors
"""

class CommandError(PotentialFittingError):
    """Exception for all called process errors to extends"""

    def __init__(self, command, message):
        super().__init__("Error when executing system command '{}': {}".format(command, message))
        self.command = command

class CommandNotFoundError(CommandError):
    """Raised when a command is not found"""

    def __init__(self, command):
        super().__init__(command, "command not found")

class CommandExecutionError(CommandError):
    """Raised when a command fails for some reason"""

    def __init__(self, command, call, exit_code, error):
        super().__init__(command, "entire command: '{}', error: {} (exit code {})".format(" ".join(call), error, exit_code))
        self.error = error
        self.call = call

    def get_error():
        return self.error

"""
--------------------------- Library Errors
"""

class LibraryError(PotentialFittingError):
    """Exception for all errors regarding the libraries we are using"""
    
    def __init__(self, library, message):
        super().__init__("Library '{}' Error: {}".format(library, message))

class LibraryCallError(LibraryError):
    """Exception for all errors caused by calls to Libraries"""

    def __init__(self, library, call, message):
        super().__init__(library, "Error during call to '{}': {}".format(call, message))

class LibraryNotAvailableError(LibraryError):
    """Exception when a user tries to use a library that is not installed"""

    def __init__(self, library):
        super().__init__(library, "library '{}' not installed or failed to be found".format(library))

"""
--------------------------- Database Errors
"""

class DatabaseError(PotentialFittingError):
    """Exception for all errors caused by database access"""
    
    def __init__(self, database, message):
        super().__init__("Error during access to database '{}': {}".format(database, message))

class InconsistentDatabaseError(DatabaseError):
    """Raised when a database contains inconsistent informaiton that should not be possible"""

    def __init__(self, database, message):
        super().__init__(database, "Inconsistent information in database: {}".format(message))

class NoEnergiesError(DatabaseError):
    """Raised when a user performs an operation that requires calculated energies when no energies are calculated"""

    def __init__(self, database, molecule_name, method, basis, cp, tag):
        super().__init__(database, "Unable to find calculated energies for molecule '{}' with the model '{}/{}' with cp={} and tag='{}'".format(molecule_name, method, basis, cp, tag))

class NoOptimizedEnergyError(DatabaseError):
    """Raised when a user performs an operation that requires an optimized energy when one does not exist"""

    def __init__(self, database, molecule_name, method, basis, cp, tag):
        super().__init__(database, "Unable to find calculated optimized energy for molecule '{}' with the model '{}/{}' with cp={} and tag='{}'".format(molecule_name, method, basis, cp, tag))

class MultipleOptimizedEnergiesError(DatabaseError):
    """Raised when a user performs an operation that requires an optimized energy when 2 or more exist"""

    def __init__(self, database, molecule_name, method, basis, cp, tag, num_opt):
        super().__init__(database, "Found multiple ({}) calculated optimized energies for molecule '{}' with the model '{}/{}' with cp={} and tag='{}'".format(num_opt, molecule_name, method, basis, cp, tag))

class NoEnergyInRangeError(DatabaseError):
    """Raised when a user performs an operation that has no valid energy configurations in the provided range of potential values."""

    def __init__(self, database, molecule_name, method, basis, cp, tag, e_min, e_max):
        super().__init__(database, "Unable to generate training set for molecule'{}' with the model '{}/{}' with e_min = '{}' and e_max = '{}'".format(molecule_name, method, basis, e_min, e_max))

"""
--------------------------- InputException
"""

class InvalidInputError(PotentialFittingError):
    """Basic exception for all invalid input from the user"""

    def __init__(self, message):
        super().__init__(message)

class XYZFormatError(InvalidInputError):
    """Raised when an XYZ file has invalid formatting"""

    def __init__(self, line, fix):
        super().__init__("Invalid xyz formatting on or near line: '{}'. Line format should be {}".format(line, fix))

class ParsingError(PotentialFittingError):
    """Exception for problems reading a file"""

    def __init__(self, file_path, message):
        super().__init__("A problem occured while parsing file '{}': {}".format(file_path, message))

class LineFormatError(ParsingError):
    """Exception for problems with the format of a line in a file"""

    def __init__(self, file_path, line, correct_format):
        super().__init__(file_path, "Error in line '{}': format should be: '{}'".format(line, correct_format))

class InvalidValueError(InvalidInputError):
    """Raised when a value a user inputs has an invalid value"""

    def __init__(self, prop, value, fix):
        super().__init__("Property '{}' has invalid value '{}'; value must be {}".format(prop, value, fix))

class InconsistentValueError(InvalidInputError):
    """Raised when input from multiple sources is inconsistent"""

    def __init__(self, property1, property2, value1, value2, fix):
        super().__init__("Properties '{}' with value '{}' and '{}' with value '{}' are inconsistent: {}".format(property1, value1, property2, value2, fix))

class NoSuchLibraryError(InvalidInputError):
    """Raised when the user requests a library that is not built in to our code"""

    def __init__(self, library):
        super().__init__("Unrecognized Library: '{}'".format(library))

class ConfigError(InvalidInputError):
    """Basic exception for all Config Errors"""

    def __init__(self, file, message):
        super().__init__("Error in settings file '{}': {}".format(file, message))

class ConfigMissingFileError(ConfigError):
    """Raise when a specified settings file does not exist"""

    def __init__(self, file):
        super().__init__(file, "file does not exist")

class ConfigMissingSectionError(ConfigError):
    """Raised when a required section is missing from a settings file"""

    def __init__(self, file, section, prop):
        super().__init__(file, "You must define section '{}' with property '{}'".format(section, prop))

class ConfigMissingPropertyError(ConfigError):
    """Raised when a required property is missing from a settings file"""

    def __init__(self, file, section, prop):
        super().__init__(file, "You must define property '{}' in section '{}'".format(prop, section))

class ConfigPropertyWrongFormatError(ConfigError):
    """Raised when a property in a settings file is in an improper format"""

    def __init__(self, section, prop, data, form):
        super().__init__("property '{}' in section '{}' is in invalid form: '{}'; correct form is: {}".format(prop, section, data, form))

class StopLoop(Exception):
    """Used as a type of break statement for nested loops"""

    def __init__(self, name):
        self.name = name
        super().__init__("Error, this exception should always be caught")
