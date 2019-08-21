#!/bin/bash

if [ $# -ne 3 ]; then
    echo "wrong number of arguments."
    echo "Usage: ./coverage.sh <test_results> <test_log> <coverage_log>"
    exit 1
fi

echo "Running python tests."

coverage run --source potential_fitting run_tests.py > $2 2> $1

coverage report -m > $3

echo "Done!" 
echo "Test results in ${1}."
echo "Test log in ${2}."
echo "Test covarege report in ${3}."
