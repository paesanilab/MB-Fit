#!/bin/bash

WORKDIR=$PWD
../run_example.sh A1B2_A1B2.in $WORKDIR/poly/ $WORKDIR/../template/ $WORKDIR/../get_2b_codes.py 

if [ $? -eq 0 ]; then
  echo "Code generation complete and compiled"
else
  echo "Code generation failed"
  exit 1
fi

cd fit
../A1B2_A1B2/fit-2b ts_co2_co2.xyz > fit.log

if [ $? -eq 0 ]; then
  echo "Fit complete"
else
  echo "Fit failed"
  exit 1
fi



