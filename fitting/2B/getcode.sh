#!/bin/bash

if [ $# -ne 4 ]; then
  echo "USAGE: $0 <arg1> <arg2> <arg3> <arg4>"
  echo "Arg 1 is the .in file used to generate the polynomials"
  echo "Arg 2 is the path to the folder were the poly*.cpp files are."
  echo "--> Assumes that Arg1 is inside the path of Arg2."
  echo "Arg 3 is the path to the files inside the template folder"
  echo "Arg 4 is the path to the python script"
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

cp $2/$1 .
cp $2/poly-direct.cpp .
cp $2/poly-grd.cpp ./poly_2b_${FNAME}_v1x.cpp
cp $2/poly-nogrd.cpp ./poly_2b_${FNAME}_v1.cpp
cp $2/poly-model.h ./poly_2b_${FNAME}_v1x.h

# CHange the includes in the poly files
cat poly_2b_${FNAME}_v1x.cpp | sed "s/poly-model.h/poly_2b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_2b_${FNAME}_v1x.cpp
cat poly_2b_${FNAME}_v1.cpp | sed "s/poly-model.h/poly_2b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_2b_${FNAME}_v1.cpp

# Execute python script
echo "Did you change the parameters for your system in the python script?"

python $4 $1 poly-direct.cpp

# Copy the rest of the files
cp $3/* .

# Compile
make clean
make

# End
