#!/bin/bash

mkdir -p gcn
mkdir -p inputs
mkdir -p logs

err=0

python3 ../src/generate_nm_configs.py settings.ini
if [ $? -ne 0 ]
then
    echo "Error: Execution of generate_nm_configs.py failed"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_optimized.xyz
diff expected/$fn inputs/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: optimized geometries differ"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_normalmodes.dat
diff expected/$fn inputs/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: normalmodes differ"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz.out
diff <(sed -e '$d' expected/$fn) <(sed -e '$d' gcn/$fn) > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: nmcgen output differs"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_configurations.xyz
diff expected/$fn gcn/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: configurations differ"
    err=$(($err+1))
fi

if [ $err -ne 0 ]
then
    echo "test failed"
else
    echo "test passed"
    rm timer.dat
    rm -r input gcn logs
    rm *.dif
fi

exit $err
