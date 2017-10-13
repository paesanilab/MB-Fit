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
#echo $nblocks
lblock=$(bc <<< "($dim/3+8)")
#echo $lblock

### record the atomic masses
##masses=$(grep -a -A $[dim/3-1] "Atom  1:" $qchem_outfile | tail -n $[dim/3] | cut -c 25-34)
##set -- $masses
##for m in `seq 1 1 $#`; do
##    mass[$m]=$1
##    shift
##    echo -e "mass[$m]: ${mass[m]}"
##done

        
# grep -a "VIBRATIONAL ANALYSIS" $qchem_outfile
k=1
j=1
for i in `seq 1 1 $nblocks`; do
    ci=12 #inital column number
    cf=31 #final column number
    until [[ ${j} -gt $nmodes ]]; do
    #for j in `seq 1 1 53`; do


        # record the normal mode coordinates reported by QChem; 
        # QChem normal mode coordinates ?are *not*? mass-scaled?
        block=$(grep -a -A $[lblock-1] "Mode:.* $k " $qchem_outfile | tail -n $[dim/3] | cut -c $ci-$cf)
        set -- $block
        for m in `seq 1 1 $dim`; do
            coord[$m]=$1
            shift
            #echo -e "coord[$m]: ${coord[m]}"
        done

        # output normal mode index, frequency, reduced mass 
        echo "normal mode:  $j"
        freq=$(grep -a -A 1 "Mode:.* $k " $qchem_outfile | tail -n 1 | cut -c $ci-$cf)
        echo $freq #frequency in cm^{-1}
        redmass=$(grep -a -A 3 "Mode:.* $k " $qchem_outfile | tail -n 1 | cut -c $ci-$cf)

        # determine the reduced mass associated with the j^{th} normal mode
##        # output the reduced mass
##        m=1; n=1
##        redmass=0
##        #for n in `seq 1 1 ${#mass[@]}`; do
##            until [[ $m -gt $dim ]]; do
##                redmass=$(bc <<< "scale=10;${coord[m]}^2/${mass[n]}+$redmass"); m=$[m+1] ##
##                redmass=$(bc <<< "scale=10;${coord[m]}^2/${mass[n]}+$redmass"); m=$[m+1] ##
##                redmass=$(bc <<< "scale=10;${coord[m]}^2/${mass[n]}+$redmass"); m=$[m+1] ##
##                n=$[n+1]
##            done
##        #done
##        redmass=$(bc <<< "scale=10;sqrt(1/$redmass)") ##
##        redmass=$(echo $redmass | awk '{printf "%12.8f\n", $1}')
        echo "red_mass =    $redmass"


##        # output the mass-scaled normal mode coordinates for the j^{th} mode; 
##        # coords[*] = coords[*]*sqrt(redmass)
##        for m in `seq 1 1 $dim`; do
##            coords[$m]=$(bc <<< "scale=10;${coord[m]}*sqrt($redmass)")
##        done

        m=1
        until [[ $m -gt $dim ]]; do
            echo ${coord[m]} | awk '{printf "   %12.3f", $1}'; m=$[m+1]
            echo ${coord[m]} | awk '{printf "   %12.3f", $1}'; m=$[m+1]
            echo ${coord[m]} | awk '{printf "   %12.3f\n", $1}'; m=$[m+1]
            #echo ${coords[m]} | awk '{printf "   %12.8f", $1}'; m=$[m+1]
            #echo ${coords[m]} | awk '{printf "   %12.8f", $1}'; m=$[m+1]
            #echo ${coords[m]} | awk '{printf "   %12.8f\n", $1}'; m=$[m+1]
        done
        printf "\n"

        ## will output the ?mass-scaled? QChem coordinates
        #for l in `seq 1 3 $dim`; do
        #    echo $block | cut -d" " -f $l-$[l+2]
        #    sleep 0.1s
        #done
        j=$[j+1]
        ci=$[ci+23]
        cf=$[cf+23]
        if [ $(( $[j-1] % 3 )) == 0 ]; then
            ci=12
            cf=31
            k=$[k+3]
        fi
    done
    #k=$[k+3]
done
