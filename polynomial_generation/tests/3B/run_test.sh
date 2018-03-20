#!/bin/bash 

generator=$(pwd)/../../src/poly-gen_mb-nrg.pl
order=2

tests="clh2o_2"

for test in $tests
do

    pushd $test > /dev/null
    $generator $order $test.in > $test.log
    diff $test.log expected/$test.log
    exit_code=$?
    if [ $exit_code -ne 0 ]
    then
	echo "$test failed!!"
	exit $exit_code
    else
	echo "$test passed!"
	rm *log poly* vars.cpp
    fi
    popd > /dev/null
    
done


