#!/bin/bash 

generator=poly-gen_mb-nrg.pl

run_test() {
    test=$1
    order=$2
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
}

order=2
test="clh2o_2"
run_test $test $order

order=4
test="mb-pol-3b"
run_test $test $order

