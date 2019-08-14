TRAVIS_BUILD_DIR=$(pwd)
num_passed=0
num_failed=0

check_exit_code () {
    if [ $? -eq 0 ]
    then
        num_passed=$(($num_passed + 1))
    else
        num_failed=$(($num_failed + 1))
    fi
}

echo "Running polynomial tests"
cd $TRAVIS_BUILD_DIR/oldtests/polytests/1B; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/polytests/2B; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/polytests/3B; bash run_test.sh
check_exit_code

echo "Running polynomial input tests"
cd $TRAVIS_BUILD_DIR/oldtests/polytests/input_generation; bash run_test.sh
check_exit_code

echo "Running calculator tests"
cd $TRAVIS_BUILD_DIR/oldtests/calctests/monomer; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/calctests/monomers; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/calctests/dimer; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/calctests/dimers; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/calctests/trimer; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/oldtests/calctests/trimers; bash run_test.sh
check_exit_code

echo "running norm distribution tests"
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/src; make
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_so4.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_pf6.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_trimer_minimum.sh
check_exit_code
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_trimer_saddle.sh
check_exit_code

echo "running psi4 tests"
cd $TRAVIS_BUILD_DIR/oldtests/calctests/psi4; bash run_test.sh
check_exit_code

cd $TRAVIS_BUILD_DIR

echo "Test running complete, $num_passed passed and $num_failed failed."
if [ $num_failed -gt 0 ]
then
    exit 1
else
    exit 0
fi