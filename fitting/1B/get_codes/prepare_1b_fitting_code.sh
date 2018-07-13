#!/bin/bash

if [ $# -ne 2 ]; then
  echo "USAGE: $0 <arg1> <arg2>"
  echo "Arg 1 is the absolute path to the .in file used to generate the polynomials"
  echo "Arg 2 is the path to the folder were the poly*.cpp files are."
  exit
fi

INPUT=$1
POLY_CPP_PATH=/server-home1/ebullvul/potential_fitting/projects/CO2/poly
PYTHON_SCRIPT=get-1b-fit.py
PYTHON_SCRIPT_PATH=/server-home1/ebullvul/potential_fitting/fitting/1B/get_codes
TEMPLATE_PATH=$PYTHON_SCRIPT_PATH/template

# Getting the system name
FNAME=$(basename "${INPUT}" .in)

cp $1 .
cp $POLY_CPP_PATH/poly-direct.cpp .
cp $POLY_CPP_PATH/poly-grd.cpp ./poly_1b_${FNAME}_v1x.cpp
cp $POLY_CPP_PATH/poly-nogrd.cpp ./poly_1b_${FNAME}_v1.cpp
cp $POLY_CPP_PATH/poly-model.h ./poly_1b_${FNAME}_v1x.h

nvar=`grep "<> variables" $POLY_CPP_PATH/poly.log | sed 's/(//g' | sed 's/)//g' | awk '{print $3}'`
npoly=`grep "Total number of terms" $POLY_CPP_PATH/poly.log | awk '{print $5}'`

cp ../config.ini ../config.initmp

cat << EOF >> ../config.ini

## Define number of variables and terms
nvars = $nvar
npoly = $npoly
EOF


# CHange the includes in the poly files
cat poly_1b_${FNAME}_v1x.cpp | sed "s/poly-model.h/poly_1b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_1b_${FNAME}_v1x.cpp
cat poly_1b_${FNAME}_v1.cpp | sed "s/poly-model.h/poly_1b_${FNAME}_v1x.h/g" > tmp
mv tmp poly_1b_${FNAME}_v1.cpp

echo "Executing python generator script"

# Execute python script
python3 $PYTHON_SCRIPT_PATH/$PYTHON_SCRIPT $INPUT poly-direct.cpp ../config.ini
mv ../config.initmp ../config.ini

# Copy the rest of the files
cp $TEMPLATE_PATH/* .

# Compile
make clean
make

# End
