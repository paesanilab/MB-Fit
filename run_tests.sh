#!/bin/bash

export MBFIT_HOME=$PWD


if [ $# -eq 0 ]
then
    time=$(date '+%Y-%m-%d-%H-%M-%S')
    mkdir -p test_results
    results="test_results/${time}_results.log"
    log="test_results/${time}_test.log"
    coverage="test_results/${time}_coverage.log"
    failed="test_results/${time}_failed_tests_outputs"
elif [ $# -eq 4 ]
then
    results=$1
    log=$2
    coverage=$3
    failed=$4
else
    echo "wrong number of arguments."
    echo "Usage: ./coverage.sh <test_results> <test_log> <coverage_log> Mfailed_test_outputs>"
    exit 1
fi

which coverage > /dev/null
if [ $? -ne 0 ]; then
    echo "coverage is not installed."
    echo "Install coverage for instance with 'conda install coverage'."
    exit 1
fi

echo "Running python tests."

coverage run --source potential_fitting run_tests.py > $log 2> $results
mv failed_tests_outputs $failed

coverage report -m > $coverage

echo "Done!" 
echo "Test results in ${results}."
echo "Test log in ${log}."
echo "Test covarege report in ${coverage}."
