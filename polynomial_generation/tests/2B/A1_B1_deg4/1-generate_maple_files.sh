#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage $0 <degree> <input.in>"
  exit
fi

PF_HOME="/home/mrierari/codes/potential_fitting/"

$PF_HOME/polynomial_generation/src/poly-gen_mb-nrg.pl $1 $2 > Deg${1}-${2}.log
