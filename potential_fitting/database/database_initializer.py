import os
import sqlite3
from .database import Database

from potential_fitting.utils import SettingsReader

from potential_fitting.exceptions import ConfigMissingSectionError, ConfigMissingPropertyError, ParsingError, InvalidValueError, XYZFormatError, InconsistentValueError

from potential_fitting.molecule import xyz_to_molecules
"""
initializes a database from the config files in a directory
"""
def initialize_database(settings_file, database_name, directory, tag = "none"):
    """
    Initializes a new database or adds energies to be calculated to an existing one.

    Args:
        settings - the .ini file used to initialize the database
        database_name - the name of the file to save the database in
        directory - the directory of the config.xyz files to add to the database or a single .xyz file

    Returns:
        None
    """

    with Database(database_name) as database:
        
        # Create the database by initializing all tables
        database.create()
    
        print("Initializing database from xyz files in {} directory into database {}".format(directory, database_name))
        # get a list of all the files in the directory, or just the filename if directory is just a single file.
        filenames = get_filenames(directory) if os.path.isdir(directory) else [directory] if directory[-4:] == ".xyz" else []

        if len(filenames) == 0:
            raise InvalidValueError("directory", directory, "a directory with 1 or more .xyz files, or a single .xyz file.")

        # parse settings.ini file
        settings = SettingsReader(settings_file)

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

        print("Initializing of database {} successful".format(database_name))


def get_filenames(directory):
    """
    gets list of all filenames in a directory

    Args:
        directory - directory to get filenames from

    Returns:
        list of all filenames in that directory, and subdirectories
    """
    filenames = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filenames.append(root + "/" + filename);
    return filenames;
