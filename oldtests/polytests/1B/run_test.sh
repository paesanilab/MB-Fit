#!/bin/bash 

generator=poly-gen_mb-nrg.pl
order=4

tests="bf4 pdh3"

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

