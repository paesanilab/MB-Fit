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

