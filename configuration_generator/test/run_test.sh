#!/bin/bash

gcn=../norm_distribution/src/generate_configs_normdistrbn

mkdir logs
mkdir diffs
mkdir diffs/logs
mkdir diffs/config_files
err=0

if [ ! -x $gcn ]
then
    "Error: You must build the norm distribution generator first"
    exit
fi

python3 ../src/nmcgen.py settings.ini h2o.xyz h2o.opt.xyz normal_modes configs.xyz
if [ $? -ne 0 ]
then
    echo "Error: Execution of generate_nm_configs.py failed"
    err=$(($err+1))
fi

diff expected/h2o.opt.xyz h2o.opt.xyz > diffs/$fn.dif
if [ $? -ne 0 ]
then
    echo "Error: optimized geometries differ - check diffs/$fn.dif"
    err=$(($err+1))
fi

diff expected/normal_modes normal_modes > diffs/$fn.dif
if [ $? -ne 0 ]
then
    echo "Error: normalmodes differ - check diffs/$fn.dif"
    err=$(($err+1))
fi

diff expected/configs.xyz configs.xyz > diffs/$fn.dif
if [ $? -ne 0 ]
then
    echo "Error: configurations differ - check diffs/$fn.dif"
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
