#!/bin/bash

if [ $# -ne 3 ]; then
  echo "USAGE: $0 <arg1> <arg2> <arg3> <arg4>"
  echo "Arg 1 is the .in file used to generate the polynomials"
  echo "Arg 2 is the path to the files inside the template folder"
  echo "Arg 3 is the path to the python script"
  exit
fi

# Getting the system name
FNAME="${1%.*}"

# Check if the system exists
if [ -d "$FNAME" ]; then
  echo "Folder $FNAME already exists."
  exit
fi

# Create folder and copy poly files to it
mkdir $FNAME
cd $FNAME

# Execute python script
echo "Did you change the parameters for your system in the python script?"

python $3 $1 

# Copy the rest of the files
cp $2/* .

# Compile
make clean
make

# End
