#!/bin/bash

driver=../../src/driver.py

python3 $driver

ndiff training_set.xyz expected/training_set.xyz > test.log
exit_code=$?

ndiff mbdecomp.log expected/mbdecomp.log >> test.log
exit_code=$(($exit_code + $?))

ndiff -quiet -abserr 1.e-12 -separators '[ \t,()]' json_output.json expected/json_output.json >> test.log
exit_code=$(($exit_code + $?))

if [ $exit_code -ne 0 ]
then
    echo "test error!"
    exit $exit_code
else
    echo "test passed!"
    rm -f timer.dat training_set.xyz mbdecomp.log json_output.json test.log
fi
