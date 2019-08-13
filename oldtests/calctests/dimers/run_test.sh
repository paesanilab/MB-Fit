#!/bin/bash

driver=../driver.py
log=test.log

python3 $driver
exit_code=$?

ndiff -quiet -abserr 1.e-5 training_set.xyz expected/training_set.xyz > $log
exit_code=$(($exit_code + $?))

ndiff -quiet -abserr 1.e-5 mbdecomp.log expected/mbdecomp.log >> $log
exit_code=$(($exit_code + $?))

ndiff -quiet -abserr 1.e-5 -separators '[ \t,()]' json_output.json expected/json_output.json >> $log
exit_code=$(($exit_code + $?))

if [ $exit_code -ne 0 ]
then
    echo "test error!"
    exit $exit_code
else
    echo "test passed!"
    rm -f timer.dat training_set.xyz mbdecomp.log json_output.json $log
    rm -rf logs
fi
