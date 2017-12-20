#!/bin/bash


../src/generate_configs_normdistrbn < ../examples/fluoride-water_trimer/saddle/input_isomer3.txt > ../examples/fluoride-water_trimer/saddle/output_isomer3.txt

# compare output files; omit the last line (which reports CPU time used) when comparing output
cat ../examples/fluoride-water_trimer/saddle/output_isomer3.txt | sed '$d' >> temp1.txt
cat ../examples/fluoride-water_trimer/saddle/expected_output/output_isomer3.txt | sed '$d' >> temp2.txt
cmp -s temp1.txt temp2.txt

exit_code1=$?

# compare configuration files
cat ../examples/fluoride-water_trimer/saddle/isomer3_linear_Q.xyz > temp1.txt
cat ../examples/fluoride-water_trimer/saddle/expected_output/isomer3_linear_Q.xyz > temp2.txt
cmp -s temp1.txt temp2.txt

exit_code2=$?


if [ $exit_code1 -ne 0 ] || [ $exit_code2 -ne 0 ]
then
    echo "Error: Output is not as expected."
    exit_code=1
else
    echo "Output is as expected."
    exit_code=0
fi

rm temp1.txt
rm temp2.txt

exit $exit_code

