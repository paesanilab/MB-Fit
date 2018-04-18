#!/bin/bash

python3 ../src/driver.py

diff training_set.xyz expected/training_set.xyz

exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "test error!"
    exit $exit_code
else
    echo "test passed!"
    rm -r timer.dat training_set.xyz
fi
