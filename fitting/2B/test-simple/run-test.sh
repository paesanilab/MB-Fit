#!/bin/bash

if [ $# -ne 3 ]; then
  echo "USAGE: $0 <input.in> <full_path_to_poly_files> <full_path_templates>"
  exit
fi

# Create new folder
filename=$1
filename="${filename##*/}"
folder="${filename%.*}"

if [ -d "$folder" ]; then
  echo "$folder already exists..."
  exit
fi

mkdir $folder
cd $folder

# Copy poly-direct
cp $2/poly-direct.cpp .

# Copy poly-grd and change the include file
tail -n +2 $2/poly-grd.cpp > tmp
echo '#include "poly_2b_'${folder}'_v1x.h"' > poly_2b_${folder}_v1x.cpp
cat tmp >> poly_2b_${folder}_v1x.cpp

# Copy poly-nogrd and change the include file
tail -n +2 $2/poly-nogrd.cpp > tmp
echo '#include "poly_2b_'${folder}'_v1x.h"' > poly_2b_${folder}_v1.cpp
cat tmp >> poly_2b_${folder}_v1.cpp

# Copy poly-model 
cp $2/poly-model.h poly_2b_${folder}_v1x.h

# Copy the templates
cp $3/* .

# Run python script
python3 ../generate_fitting_notebook.py $1 poly-direct.cpp
rm poly-direct.cpp

# Compile
make
