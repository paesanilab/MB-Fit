#!/bin/bash

python3 ../scripts/run_psi4_1B.py

diff training_set.xyz expected/training_set.xyz

exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "run_psi4_1b.py error!"
    exit $exit_code
else
    echo "run_psi4_1b.py passed!"
    rm -r calculations* timer.dat training_set.xyz
fi
