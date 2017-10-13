#!/bin/bash


../src/generate_configs_normdistrbn < ../examples/hexafluorophosphate/input_pf6.txt > ../examples/hexafluorophosphate/output_pf6.txt

# omit the last line (which reports CPU time used) when comparing output
cat ../examples/hexafluorophosphate/output_pf6.txt | sed '$d' >> temp1.txt
cat ../examples/hexafluorophosphate/expected_output/output_pf6.txt | sed '$d' >> temp2.txt
cmp -s temp1.txt temp2.txt

exit_code=$?

if [ $exit_code -ne 0 ]
then
    echo "Error: Output is not as expected."
else
    echo "Output is as expected."
fi

rm temp1.txt
rm temp2.txt

exit $exit_code

