#!/bin/bash
for j in $1/*
do
    if [[ -f $j ]]
    then
        python3 $j
    fi
done