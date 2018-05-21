"""
DEFINES EXCEPTIONS NEEDED BY THE PROGRAM
"""

"""
Raised when the user provides input file in incorrect format
"""
class InvalidFormatException(Exception):
    def __init__(self, f, line, message):
        self.message = "Invalid formatting of input file {} at line {}: {}".format(f, line, message)

"""
Raised when an Atom's property is set to an ivalid number
"""
class InvalidPropertyException(Exception):
    def __init__(self, message):
        self.message = message

"""
Raised when the user does not provide a database
"""
class DatabaseUnspecifiedException(Exception):
    def __init__(self, config, location):
        self.message = "Must specify name of database in {} under {}".format(config, location)

"""
Raised when a database operation is performed without first setting the id
"""
class MoleculeIDUnsetException(Exception):
    def __init__(self, operation, database):
        self.message = "Must first set molecule ID with Database.set_molecule_id(SHA1_id) before performing operation {} on database {}".format(operation, database)

"""
Raised when the user tries to modify a row in a database without first initializing it
"""
class RowOperationBeforeInitException(Exception):
    def __init__(self, operation, database):
        self.message = "Must class init_molecule() before performing operation {} in database {}".format(operation, database)
