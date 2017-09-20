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

INPUT=$1
POLY_CPP_PATH=$2
TEMPLATE_PATH=$3
PYTHON_SCRIPT=$4

# Getting the system name
FNAME="${INPUT%.*}"

# Check if the system exists
if [ -d "$FNAME" ]; then
  echo "Folder $FNAME already exists."
  exit
fi

# Create folder and copy poly files to it
mkdir $FNAME
cd $FNAME

cp $POLY_CPP_PATH/$1 .
cp $POLY_CPP_PATH/poly-direct.cpp .
cp $POLY_CPP_PATH/poly-grd.cpp ./poly_1b_${FNAME}_v1x.cpp
cp $POLY_CPP_PATH/poly-nogrd.cpp ./poly_1b_${FNAME}_v1.cpp
cp $POLY_CPP_PATH/poly-model.h ./poly_1b_${FNAME}_v1x.h

# CHange the includes in the poly files
cat poly_1b_${FNAME}_v1x.cpp | sed "s/poly-model.h/poly_1b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_1b_${FNAME}_v1x.cpp
cat poly_1b_${FNAME}_v1.cpp | sed "s/poly-model.h/poly_1b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_1b_${FNAME}_v1.cpp

# Execute python script
echo "Did you change the parameters for your system in the python script?"

python $PYTHON_SCRIPT $INPUT poly-direct.cpp

# Copy the rest of the files
cp $TEMPLATE_PATH/* .

# Compile
make clean
make

# End
