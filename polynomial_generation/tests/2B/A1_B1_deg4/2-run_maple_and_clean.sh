#!/bin/bash

PF_HOME="/home/mrierari/codes/potential_fitting/"

maple poly-grd.maple
maple poly-nogrd.maple

$PF_HOME/polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp
$PF_HOME/polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp

