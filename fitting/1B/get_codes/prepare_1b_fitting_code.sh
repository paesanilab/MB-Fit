#!/bin/bash

if [ $# -ne 4 ]; then
  echo "USAGE: $0 <arg1> <arg2> <arg3> <arg4>"
  echo "Arg 1 is the absolute path to the .in file used to generate the polynomials"
  echo "Arg 2 is the path to the folder were the poly*.cpp files are."
  echo "Arg 3 is the path to the files inside the template folder"
  exit
fi

INPUT=$1
POLY_CPP_PATH=../polynomial_generation
TEMPLATE_PATH=$3
PYTHON_SCRIPT=get-1b-fit.py

# Getting the system name
FNAME=$(basename "${INPUT}")

cp $1 .
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

$PYTHON_SCRIPT $INPUT poly-direct.cpp ../config.ini

# Copy the rest of the files
cp $TEMPLATE_PATH/* .

# Compile
make clean
make

# End
