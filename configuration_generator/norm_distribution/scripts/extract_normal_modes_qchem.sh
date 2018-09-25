#!/bin/bash

# QChem normal modes do not have mass-scaling included
# QChem atomic coordinates given in bohr by default

if [[ $# -ne 3 ]]; then
    echo "Usage: sh extract_normalmodes.sh <dim> <dimnull> <QChem_outfile>"
    exit 1
fi

dim=$1
dimnull=$2
qchem_outfile=$3

nmodes=$(bc <<< "($dim - $dimnull)/1")
nblocks=$(bc -l <<< "($nmodes/3)")
nblocks=$(bc <<< "($nblocks+0.5)/1")
lblock=$(bc <<< "($dim/3+8)")

        
# grep -a "VIBRATIONAL ANALYSIS" $qchem_outfile
k=1
j=1
for i in `seq 1 1 $nblocks`; do
    ci=12 #inital column number
    cf=31 #final column number
    until [[ ${j} -gt $nmodes ]]; do

        # record the normal mode coordinates reported by QChem:
        block=$(grep -a -A $[lblock-1] "Mode:.* $k " $qchem_outfile | tail -n $[dim/3] | cut -c $ci-$cf)
        set -- $block
        for m in `seq 1 1 $dim`; do
            coord[$m]=$1
            shift
        done

        # output normal mode index, frequency, reduced mass 
        echo "normal mode:  $j"
        freq=$(grep -a -A 1 "Mode:.* $k " $qchem_outfile | tail -n 1 | cut -c $ci-$cf)
        echo $freq #frequency in cm^{-1}
        redmass=$(grep -a -A 3 "Mode:.* $k " $qchem_outfile | tail -n 1 | cut -c $ci-$cf)
        echo "red_mass =    $redmass"

        # output normal mode coordinates
        m=1
        until [[ $m -gt $dim ]]; do
            echo ${coord[m]} | awk '{printf "   %12.3f", $1}'; m=$[m+1]
            echo ${coord[m]} | awk '{printf "   %12.3f", $1}'; m=$[m+1]
            echo ${coord[m]} | awk '{printf "   %12.3f\n", $1}'; m=$[m+1]
        done
        printf "\n"

        j=$[j+1]
        ci=$[ci+23]
        cf=$[cf+23]
        if [ $(( $[j-1] % 3 )) == 0 ]; then
            ci=12
            cf=31
            k=$[k+3]
        fi
    done
done
