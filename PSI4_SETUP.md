# Psi4 Setup Guide

# General:
This guide will show you how to setup your anaconda enviornment to run the psi4 development branch.

# Install anaconda

See https://conda.io/docs/user-guide/install to install anaconda. If given the choice, install full Anaconda, not Miniconda

Update anaconda to make sure everything installed properly:

```
conda update conda
```

# Create an enviornment

Create a new conda enviornment:
```
conda create -n psi4dev python=3.6 anaconda
```

The argument following -n specifies the name of the enviornment, feel free to chose whatever you like.

python=#.# specifies the version of python to use. Currently, psi4 works best with python 3.6

# Activate your enviornment

Activate your new enviornment:
```
conda activate psi4dev
```

If you used a name other than psi4dev when creating the enviornment, you will have to use that name here in place of psi4dev,

# Install the psi4 development branch

Install the development branch of psi4 in your anaconda enviornment:
```
conda install -c psi4/label/dev psi4
```

# You're done!

Your acaconda enviornment should be set up to run the psi4 development branch.

Remember, you can deactivate your enviornment with
```
conda deactivate
```
And you can activate your enviornment with
```
conda activate psi4dev
```
