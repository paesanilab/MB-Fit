from database import Database
from molecule_parser import xyz_to_molecules
import configparser
import os

try:
    os.unlink("testdb.db")
except FileNotFoundError:
    pass

# init a molecule
f = open("config_files/trimer.xyz", 'r')

# get list of all files and subdirectories in directory
config = configparser.SafeConfigParser(allow_no_value=False)
# read the config file into the configuration
config.read("config_files/settings.ini")
molecule = xyz_to_molecules(f, config)[0]
print("Trimer Test:")

# open the database
database = Database("testdb.db")

# set the database's Id
database.set_molecule_id(molecule.get_SHA1())

# initialize database table rows
database.init_molecule("CONFIG_TEMP", molecule.get_num_atoms(), molecule.get_num_fragments(), "sometag", "somemodel", "cp boolean")

# check if an energy value is in the table
print("E12 is in table: " + str(database.has_nmer_energy([1, 2])))

# add an energy to the table
database.set_nmer_energy([1, 2], 5)

# check if an energy value is in the table
print("E12 is in table: " + str(database.has_nmer_energy([1, 2])))

# check if another energy value is in the table
print("E0 is in table: " + str(database.has_nmer_energy([0])))

# add another energy to the table
database.set_nmer_energy([0], -12)

# retrieve energy from table
print("E0 energy: " + str(database.get_nmer_energy([0])))

# retrieve another energy from table
print("E12 energy: " + str(database.get_nmer_energy([1, 2])))

# retrieve another energy from table
print("E012 energy: " + str(database.get_nmer_energy([0, 1, 2])))

# init a new molecule
f = open("config_files/dimers.xyz", "r")

# get list of all files and subdirectories in directory
config = configparser.SafeConfigParser(allow_no_value=False)
# read the config file into the configuration
config.read("config_files/settings_dimer.ini")
molecule = xyz_to_molecules(f, config)[0]
print("Dimer Test:")
# reset the database's Id
database.set_molecule_id(molecule.get_SHA1())

# initialize database table rows
database.init_molecule("CONFIG_TEMP", molecule.get_num_atoms(), molecule.get_num_fragments(), "tag2", "modelmethod", "cp boolean")

# check if an energy is in the table
print("E12 is in table: " + str(database.has_nmer_energy([1, 2])))

# add an energy to the table
database.set_nmer_energy([1, 2], 5)

# check if an energy value is in the table
print("E12 is in table: " + str(database.has_nmer_energy([1, 2])))

# check if another energy value is in the table
print("E0 is in table: " + str(database.has_nmer_energy([0])))

# add another energy to the table
database.set_nmer_energy([0], -12)

# retrieve energy from table
print("E0 energy: " + str(database.get_nmer_energy([0])))

# retrieve another energy from table
print("E12 energy: " + str(database.get_nmer_energy([1, 2])))

# retrieve another energy from table
print("E012 energy: " + str(database.get_nmer_energy([0, 1, 2])))


# close the database
database.close()
