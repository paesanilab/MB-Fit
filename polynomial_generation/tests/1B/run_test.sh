#!/bin/bash 

./../../src/poly-gen_mb-nrg.pl 4 bf4.in > bf4.log 
## Testing BF4
diff bf4.log bf4_expected/bf4.log  
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "bf4 failed!!"
    exit $exit_code
else
    echo "bf4 passed!"
    rm *log poly* vars.cpp
fi

./../../src/poly-gen_mb-nrg.pl 4 pdh3.in > pdh3.log 
## Testing PdH3
diff pdh3.log pdh3_expected/pdh3.log  
exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "pdh3 failed!!"
    exit $exit_code
else
    echo "pdh3 passed!"
    rm *log poly* vars.cpp
fi

