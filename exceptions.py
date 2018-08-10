"""
--------------------------- General Exceptions
"""

class PotentialFittingError(Exception):
    """Basic exception for all errors raised by our code"""
    
    def __init__(self, message):
        super().__init__("The Following Error Occured in the Potential Fitting Library: {}".format(message))

"""
--------------------------- Library Errors
"""

class LibraryError(PotentialFittingError):
    """Exception for all errors regarding the libraries we are using"""
    
    def __init__(self, library, message):
        super().__init__("Library {} Error: {}".format(library, message))

class LibraryCallError(LibraryError):
    """Exception for all errors caused by calls to Libraries"""

    def __init__(self, library, call, message):
        super().__init__(library, "Error during call to {}: {}".format(call, message))

class LibraryNotAvailableError(LibraryError):
    """Exception when a user tries to use a library that is not installed"""

    def __init__(self, library):
        super().__init__(library, "library {} not installed or failed to be found".format(library))

"""
--------------------------- Database Errors
"""

class DatabaseError(PotentialFittingError):
    """Exception for all errors caused by database access"""
    
    def __init__(self, database, message):
        super().__init__("Error during access to database {}: {}".format(database, message))

class InconsistentDatabaseError(DatabaseError):
    """Raised when a database contains inconsistent informaiton that should not be possible"""

    def __init__(self, database, message):
        super().__init__(database, "Inconsistent information in database: {}".format(message))

"""
--------------------------- InputException
"""

class InvalidInputError(PotentialFittingError):
    """Basic exception for all invalid input from the user"""

    def __init__(self, message):
        super().__init__(message)

class XYZFormatError(InvalidInputError):
    """Raised when an XYZ file has invalid formatting"""

    def __init__(self, path, message):
        super().__init__("Invalid xyz file formatting in file {}: {}".format(path, message))

class NoSuchLibraryError(InvalidInputError):
    """Raised when the user requests a library that is not built in to our code"""

    def __init__(self, library):
        super().__init__("Unrecognized Library: {}".format(library))

class ConfigError(InvalidInputError):
    """Basic exception for all Config Errors"""

    def __init__(self, message):
        super().__init__("Error in settings file: {}".format(message))

class ConfigMissingSectionError(ConfigError):
    """Raised when a required section is missing from a settings file"""

    def __init__(self, section, prop):
        super().__init__("You must define section {} with property {}".format(section, prop))

class ConfigMissingPropertyError(ConfigError):
    """Raised when a required property is missing from a settings file"""

    def __init__(self, section, prop):
        super().__init__("You must define property {} in section {}".format(prop, section))

class StopLoop(Exception):
    """Used as a type of break statement for nested loops"""

    def __init__(self, name):
        self.name = name
        super().__init__("Error, this exception should always be caught")
