#!/bin/bash

WORKDIR=$PWD

cd A1B4
rm -rf A1B4
./run_example.sh A1B4.in $WORKDIR/A1B4/A1B4_poly $WORKDIR/template $WORKDIR//get-1b-fit.py

cd $WORKDIR/A1B4/A1B4_fit
rm fit-1b.nc 
ncgen -o fit-1b.nc fit-1b.cdl
../A1B4/eval-1b fit-1b.nc test.xyz > test_eval.out
if diff test_eval.out test_eval.out_expected &> /dev/null ; then
    echo "Test for PASSED"
  else
    echo "Test for FAILED"
    echo "EXPECTED:"
    cat test_eval.out_expected
    echo
    echo "OUTPUT:"
    cat test_eval.out
    exit 1
fi


