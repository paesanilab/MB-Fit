## POSTGRES DATBASE NOTES:

For now, please ALWAYS include your name or initials in the tag. This is so your training and test sets don't get mixed up with other people's.

# MB-nrg potential energy function generation

[![Build Status](https://travis-ci.org/paesanilab/potential_fitting.svg?branch=master)](https://travis-ci.org/paesanilab/potential_fitting)

This repository contains all of our codes to generate MB-nrg potential
energy functions (1B, 2B, 3B):

1. Training set generator consisting of
   - Configuration generator
   - QM MB energy calculator
2. Potential energy function generator consisting of
   - Code to generate polynomials functions
   - Fitting codes


## Single script

`generate_1b_fitcode.py` is a single script that executes all the steps of polynomial fitting.

First create a folder, as an example see the `B1F4` folder in the repository.
Inside this folder, create a configuration file `config.ini` that contains all
the configuration options. See for example `B1F4/config.ini`.

Then we assume that all of the tools in this repository will be installed and
available in the path. Currently we do not install them so we modify the path
by running:

    source set_path.sh

Finally we can run the script:

    cd B1F4
    python ../generate_1b_fitcode.py config.ini
    
If there are problems with the libraries:
1. Make sure you are running python3
2. Remove all python path libraries by doing `export PYTHONPATH=`
