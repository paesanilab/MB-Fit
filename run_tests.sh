TRAVIS_BUILD_DIR=$(pwd)

echo "Running polynomial tests"
cd $TRAVIS_BUILD_DIR/oldtests/polytests/1B; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/polytests/2B; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/polytests/3B; bash run_test.sh

echo "Running polynomial input tests"
cd $TRAVIS_BUILD_DIR/oldtests/polytests/input_generation; bash run_test.sh

echo "Running calculator tests"
cd $TRAVIS_BUILD_DIR/oldtests/calctests/monomer; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/calctests/monomers; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/calctests/dimer; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/calctests/dimers; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/calctests/trimer; bash run_test.sh
cd $TRAVIS_BUILD_DIR/oldtests/calctests/trimers; bash run_test.sh

echo "running norm distribution tests"
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/src; make
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test.sh
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_so4.sh
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_pf6.sh
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_trimer_minimum.sh
cd $TRAVIS_BUILD_DIR/configuration_generator/norm_distribution/scripts; bash run_test_trimer_saddle.sh

echo "running psi4 tests"
cd $TRAVIS_BUILD_DIR/oldtests/calctests/psi4; bash run_test.sh
