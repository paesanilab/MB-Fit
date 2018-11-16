# external package imports
import os, sqlite3

# absolute module imports
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import (ConfigMissingSectionError, ConfigMissingPropertyError, ParsingError,
        InvalidValueError, XYZFormatError, InconsistentValueError)
from potential_fitting.molecule import xyz_to_molecules

# local module imports
from .database import Database

def initialize_database(settings_path, database_path, directory_path, tag = "none"):
    """
    Initializes a new database or adds energies to be calculated to an existing one.

    Args:
        settings_path       - Local path to the ".ini" file with all relevant settings.
        database_path       - Local path to the database file. ".db" will be appended if it does not already end in
                ".db".
        directory_path      - Local path to the directory of ".xyz" files or a single ".xyz" file to add to the
                database. If a directory is specified, all non ".xyz" or ".xyz.opt" files will be ignored"

    Returns:
        None.
    """

    with Database(database_path) as database:
        
        # Create the database by initializing all tables
        database.create()
    
        print("Initializing database from xyz files in {} directory into database {}".format(directory_path, database_path))
        # get a list of all the files in the directory, or just the filename if directory is just a single file.
        filenames = (get_filenames(directory_path) if os.path.isdir(directory_path) else
                [directory_path] if directory_path[-4:] == ".xyz" else [])

        if len(filenames) == 0:
            raise InvalidValueError("directory", directory, 
                    "a directory with 1 or more .xyz files, or a single .xyz file.")

        # parse settings.ini file
        settings = SettingsReader(settings_path)

        # Read values from the settings file
        method = settings.get("energy_calculator", "method")

        basis = settings.get("energy_calculator", "basis")

        cp = settings.getboolean("energy_calculator", "cp")

        # loop thru all files in directory
        for filename in filenames:

            # if the filename does not end in .xyz, then skip it
            if filename[-4:] != ".xyz":
                continue

            # get list of all molecules in file
            try:
                molecules = xyz_to_molecules(filename, settings)
            except (XYZFormatError, InconsistentValueError) as e:
                raise ParsingError(filename, str(e)) from e

            # does this file contain optimized geometries?
            optimized = filename.endswith(".opt.xyz")
            
            # for each molecule in the file
            for molecule in molecules:

                molecule.move_to_center_of_mass()
                molecule.rotate_on_principal_axes()

                # add this molecule to the database
                database.add_calculation(molecule, method, basis, cp, tag, optimized)

        print("Initializing of database {} successful".format(database_path))


def get_filenames(directory):
    """
    Gets list of all filenames in a directory and any subdirectories.

    Args:
        directory           - Directory to get filenames from.

    Returns:
        List of all filenames in that directory, and any subdirectories.
    """
    filenames = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filenames.append(root + "/" + filename);
    return filenames;
