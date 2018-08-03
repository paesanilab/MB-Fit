#!/bin/bash

# Code: ../../bin/2b_ts_rigid
# Argument 1: MD trajectory of a single molecule of monomer 1
# The center of mass of this one will be kept at 
# the origin of coordinates (0.0,0.0,0.0)

# Argument 2: MD trajectory of a single molecule of monomer 2
# This one will be moved away from mol1

# Argument 3: Number of configurations per distance

# Argument 4: Minimum acceptable distance between any 
# pair of atoms AB, where A is from mol1 and B is from mol2

# Argument 5: Random number seed

#######################
# NOTES
#######################
# The step and maximum distances are hardcoded in the source.
# Usually, you should not change them, but might be of interest to have 
# only configurations in the short range or mid range.

co21_4000K.xyz  run.sh  wat1_1200K.xyz

../../bin/2b_ts_flex wat1_1200K.xyz co21_4000K.xyz 2 1.2 8436
