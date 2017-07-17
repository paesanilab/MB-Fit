#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage $0 TS"
  exit
fi

nl=`cat $1 | wc -l`
nf=$(($nl / 8))

rm -f ts_reordered.xyz
for i in `seq 1 1 $nf`; do
  head -n $(($i * 8)) $1 | tail -n 8 > tmp
  head -n 2 tmp >> ts_reordered.xyz
  tail -n 3 tmp >> ts_reordered.xyz
  tail -n 6 tmp | head -n 3 >> ts_reordered.xyz
  rm tmp
done
