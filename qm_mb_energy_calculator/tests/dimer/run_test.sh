#!/bin/bash
export PATH=$PATH:../../../tools/ndiff-2.00

driver=../../src/driver.py
log=test.log

python3 $driver

echo "training_set.xyz\n" > $log

ndiff training_set.xyz expected/training_set.xyz >> $log
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "training_set.xyz differs!"
fi

echo "mbdecomp.log\n" >> $log

ndiff mbdecomp.log expected/mbdecomp.log >> $log
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "mbdecomp.log differs!"
fi

echo "json_output.json\n" >> $log

ndiff -quiet -abserr 1.e-8 -separators '[ \t,()]' json_output.json expected/json_output.json >> $log
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "json_output.json differs!"
fi
