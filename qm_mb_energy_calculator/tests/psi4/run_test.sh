#!/bin/bash
export PATH=$PATH:../../../tools/ndiff-2.0

out=psi4_test.out

python3 psi4_test.py > $out

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
