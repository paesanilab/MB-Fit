#!/bin/bash
export PATH=$PATH:../../../tools/ndiff-2.00

out=sand_test.out

python3 generate_configs_normdistrbn.py<./input.txt> > $out

ndiff -quiet -abserr 1.e-8 $out expected/$out
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "test error!"
    exit $exit_code
else
    echo "test passed!"
    rm -f timer.dat psi4.out $out
fi
