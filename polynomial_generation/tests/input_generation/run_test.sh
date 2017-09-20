#!/bin/bash

a="A1B2_A1B2.in"
b="A1B2Z2_D1E2.in"
c="A1B4.in"

EXIT_CODE=0
for i in $a $b $c ; do
  python3 ../../src/generate_input_poly.py $i
  if diff ${i}_expected $i &> /dev/null ; then
    echo "Test for $i PASSED"
  else
    echo "Test for $i FAILED"
    echo "EXPECTED:"
    cat ${i}_expected
    echo
    echo "OUTPUT:"
    cat $i 
    EXIT_CODE=1
  fi
done

exit $EXIT_CODE
