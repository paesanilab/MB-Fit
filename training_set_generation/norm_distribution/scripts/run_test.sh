#!/bin/bash


../src/generate_configs_normdistrbn < ../examples/sulfate/input_so4.txt > ../examples/sulfate/output_so4.txt

# compare output files; omit the last line (which reports CPU time used) when comparing output
cat ../examples/sulfate/output_so4.txt | sed '$d' >> temp1.txt
cat ../examples/sulfate/expected_output/output_so4.txt | sed '$d' >> temp2.txt
cmp -s temp1.txt temp2.txt

exit_code1=$?


# compare configuration files; omit the last line (which reports CPU time used) when comparing output
cat ../examples/sulfate/so4_linear_Q.xyz | sed '$d' > temp1.txt
cat ../examples/sulfate/expected_output/so4_linear_Q.xyz | sed '$d' > temp2.txt
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

