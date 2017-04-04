#!/bin/bash 

./../../src/poly-gen_mb-nrg.pl 2 clh2o_2.in > clh2o_2.log 
## Testing Cl-(H2O)2
diff clh2o_2.log clh2o_2_expected/clh2o_2.log  

if [ $? -ne 0 ]
then
    echo "Cl-(H2O)2 failed!!"
    exit
else
    echo "Cl-(H2O)2 passed!"
    rm *log poly* vars.cpp
fi

