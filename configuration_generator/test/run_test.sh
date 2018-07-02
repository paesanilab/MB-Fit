#!/bin/bash

gcn=../norm_distribution/src/generate_configs_normdistrbn

mkdir -p outputs
err=0

if [ ! -x $gcn ]
then
    "Error: You must build the norm distribution generator first"
    exit
fi

python3 ../src/nmcgen.py settings.ini
if [ $? -ne 0 ]
then
    echo "Error: Execution of generate_nm_configs.py failed"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_optimized.xyz
diff expected/$fn outputs/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: optimized geometries differ - check $fn.dif"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_normalmodes.dat
diff expected/$fn outputs/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: normalmodes differ - check $fn.dif"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_gcn.out
diff <(sed -e '$d' expected/$fn) <(sed -e '$d' outputs/$fn) > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: nmcgen output differs - check $fn.dif"
    err=$(($err+1))
fi

fn=h2o_psi4_blyp_cc-pvdz_configurations.xyz
diff expected/$fn outputs/$fn > $fn.dif
if [ $? -ne 0 ]
then
    echo "Error: configurations differ - check $fn.dif"
    err=$(($err+1))
fi

if [ $err -ne 0 ]
then
    echo "test failed"
else
    echo "test passed"
    ./clean_test.sh
fi

exit $err
