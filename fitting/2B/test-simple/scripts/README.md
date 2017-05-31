This folder contains the python script in charge of generating the system dependant files, and the bash script that will create the folder where the fitting code will be put, and compile it.

The python script is the same as the one in the notebook. This generates the system dependent files, such as the different monomer classes, the structures taht will call the polynomials, and some other relevant cpp/h files.

The bash script copies the common files from the template folder, runs the python script, and compiles the code.
