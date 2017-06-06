#!/bin/bash

if [ $# -ne 3 ]; then
  echo "USAGE: $0 <input.in> <full_path_to_poly_files> <full_path_templates>"
  exit
fi

# Create new folder
# Gets the folder name from the input 
# in /home/mrierari/codes/potential_fitting/polynomial_generation/tests/2B/A1B2_A1B2/A1B2_A1B2.in
filename_wpath=$1
path_to_poly_files=$2
path_to_templates=$3


filename="${filename_wpath##*/}"
folder="${filename%.*}"

if [ -d "$folder" ]; then
  echo "$folder already exists..."
  exit
fi

mkdir $folder
cd $folder

# Copy poly-direct
cp $path_to_poly_files/poly-direct.cpp .

# Copy poly-grd and change the include file
tail -n +2 $path_to_poly_files/poly-grd.cpp > tmp
echo '#include "poly_2b_'${folder}'_v1x.h"' > poly_2b_${folder}_v1x.cpp
cat tmp >> poly_2b_${folder}_v1x.cpp

# Copy poly-nogrd and change the include file
tail -n +2 $path_to_poly_files/poly-nogrd.cpp > tmp
echo '#include "poly_2b_'${folder}'_v1x.h"' > poly_2b_${folder}_v1.cpp
cat tmp >> poly_2b_${folder}_v1.cpp

# Copy poly-model 
cp $path_to_poly_files/poly-model.h poly_2b_${folder}_v1x.h

# Copy the templates
cp $path_to_templates/* .

# Run python script
python3 ../scripts/generate_fitting_notebook.py $filename_wpath poly-direct.cpp
rm poly-direct.cpp

# Compile
make
