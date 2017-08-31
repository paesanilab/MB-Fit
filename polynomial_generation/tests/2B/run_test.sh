#!/bin/bash 

./../../src/poly-gen_mb-nrg.pl 2 h2o_2.in > h2o_2.log 
## Testing (H2O)2
diff h2o_2.log h2o_2_expected/h2o_2.log  
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "(H2O)2 failed!!"
    exit $exit_code
else
    echo "(H2O)2 passed!"
    rm *log poly* vars.cpp
fi

./../../src/poly-gen_mb-nrg.pl 2 n2o5_h2o.in > n2o5_h2o.log 
## Testing N2O5-H2O
diff n2o5_h2o.log n2o5_h2o_expected/n2o5_h2o.log  
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "N2O5-H2O failed!!"
    exit $exit_code
else
    echo "N2O5-H2O passed!"
    rm *log poly* vars.cpp
fi

