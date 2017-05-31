#!/bin/bash

# Getting full paths of the relative paths, since the run uses the full ones...
cd ../../../polynomial_generation/tests/2B/A1B2_A1B2/
path1=$PWD
cd -

cd ../../../polynomial_generation/tests/2B/A1B2_A1B2/
path2=$PWD
cd -

cd ../template/
path3=$PWD
cd -

# Linking the python script here, since the run.sh script 
# assumes the script is in the same folder where the test is run
ln -s scripts/generate_fitting_notebook.py .

# Running the fitting code generator script
./scripts/run.sh $path1/A1B2_A1B2.in $path2 $path3

# Remove the link to clean up things
rm generate_fitting_notebook.py

# Running the fitting code
cd CO2-CO2
echo "Running fit..."
../A1B2_A1B2/fit-2b training_set_mbnrg.xyz > fit.log

######
## Note: running the code twice will yield different results, 
## since the initial guess is randomized
######



