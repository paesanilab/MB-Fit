# Training Set Generation

## Overview:
This program can be used to calculate the many-body decomposition energies, as well as the k-mer energies of a molecule.
These molecules can also be defined as fragments by the user, inside an xyz input file.
As of June 2018, we are attempting to integrate a database into the program to streamline calcuations.

## To use this program:

### Setup:
To begin, it is assumed that the user has installed a `conda` build of some kind.
Otherwise, please follow this link: https://conda.io/docs/user-guide/install/index.html

To start, create a `conda` environment that includes Python 3, and TensorFlow. If TensorFlow is already installed, then update TensorFlow.
```
conda create -n foo python=3.6 tensorflow
```

The above command creates an environment named `foo` with Python 3.6, as well as TensorFlow 1.4.0. If your system has Python 3 by default, then feel free to exclude the Python argument like so:
```
conda create -n foo tensorflow
```

To start, create a `conda` environment that includes Python 3. The following command creates a conda environment named `foo` that has Python 3.6:
```
conda create -n foo python=3.6
```

Once you have created your environment, activate it with this command:
```
source activate foo
```

At this point, you should be in the environment `foo`.

Next, install `psi4`. The command is:
```
conda install -c psi4 psi4
```

Once you have finished installing `psi4`, feel free to run the test script under the tests directory:
```
cd tests
./test.sh
```

### To run this program:

For this program, you can run the program anywhere, as long as you have a `settings.ini` file and an xyz input file in your current working directory. The `input` option in `settings.ini` must match the name of the xyz file, including file extension.

The xyz input follows standard xyz file conventions, but we require the comment line to include a list a numbers to denote the number of atoms per fragment. For example, for a water trimer:
```
9
3 3 3
O        -0.0005900753    0.0000000000   -0.0004526740
H         0.9256102457    0.0000000000   -0.2481132352
H        -0.0000930743    0.0000000000    0.9582869092
O        -0.0005900753    2.0000000000   -0.0004526740
H         0.9256102457    2.0000000000   -0.2481132352
H        -0.0000930743    2.0000000000    0.9582869092
O        -0.0005900753    4.0000000000   -0.0004526740
H         0.9256102457    4.0000000000   -0.2481132352
H        -0.0000930743    4.0000000000    0.9582869092
```

To distinguish terminal commands and output, a `$` is being used here. To run this program:
```
$ ls
settings.ini 2mol_trimer.xyz
$ python ../src/driver.py
```

## Database usage:
We assume that you have created a `settings.ini` with the specified settings, as well as a `.xyz` input file. At the bare minimum, we ask that all values in [molecule] to be defined.
First, run `python database_initializer.py` on the directory where your `settings.ini` is, which should be your working directory.
Once completed, run `python database_filler.py` on the same directory. You may need `psi4` when running this script.

## About `settings.ini`:

Use this file to customize your options to use this program. This file must be in the same directory as your current working directory.

#### [driver]: Settings pertaining to the initialization of the program.
* input: The name of the xyz input file, including its file extension.
* model: The calculation package we are using for calculations. Currently only supports `psi4`.

#### [molecule]: Information about the molecule. This must be manually defined by the user.
* fragments: A comma-separated list to indicate the number of atoms in each fragment. By default, we assume the configuation is a monomer.
* charges: A comma-separated list to indicate the charge of each fragment. By default, all fragments have neutral charge.
* spin: A comma-separated list to indicate the spin multiplicity of each fragment. By default, all fragments have a spin multiplicity of 1.

#### [model]: Information about the methods regarding calculation.
* method: The method used for calculation.
* basis: The basis set used for calculation.
* cp: A boolean to whether use counterpoise correction in calculations.

#### [psi4]: Settings specific to the `psi4` model.
* memory: The memory `psi4` will use for calculation.
* threads: The amount of parallel threads `psi4` will use for calculation.


#### [MBdecomp]: Settings pertaining many-body analysis of the molecule.
* mbdecomp: Flag for enabling many-body decomposition.
* cluster_analysis: Flag for writing output for cluster analysis. This must be set to `True` for the bottom few flags to work.
* overwrite: Determines whether we wish to overwrite an output file, if one already exists. If append is also set to true, then overwrite will not occur.
* append: Determines whether we wish to append to an output file, if one already exists.
* max_nbody_energy: Prints up to the N-body energy of the molecule specified by this setting.
* kbody_energy: Flag for enabling the many-body decomposition of the fragment subsets of the molecule.

## To use TensorMol:
Simply follow the instructions given on their repository here: https://github.com/jparkhill/TensorMol

When performing `pip install`, use the command `pip install -e .` instead (no `sudo` permissions needed!). Once finished, download their pre-trained neural network to perform tests.
