#### Training Set Generation
This directory contains the scripts to read in geometries using Open Babel's `pybel` and calculates their 1B energies using [PSI4](http://www.psicode.org), an open-source suite of ab initio quantum chemistry programs designed for efficient, high-accuracy simulations of a variety of molecular properties.

Open Babel's `openbabel` module provides direct access to the C++ Open Babel library from Python [1]. The module comes with the `pybel` wrapper which provides convenience functions and classes that make it simpler to use the Open Babel libraries from Python, especially for file input/output and for accessing the attributes of atoms and molecules [2].  

Inside the `scripts` directory, `run_psi4_1B.py` uses `pybel` to read in a one-body xyz file and generate a one-body training set using PSI4. PSI4 output files are placed in a `calculations` directory within the current directory and are compressed into a zip file.
 
In the future, this will contain a script(s) to read in geometries and calculate 1B, 2B, 3B, and NB reference energies using PSI4. 

References:

1.[Noel M. O’Boyle, Michael Banck, Craig A. James, Chris Morley, Tim Vandermeersch, Geoffrey R. Hutchison. Open Babel: An open chemical toolbox. J. Cheminf. 2011, 3, 33.](https://jcheminf.springeropen.com/articles/10.1186/1758-2946-3-33)

2.[N.M. O’Boyle, C. Morley and G.R. Hutchison. Pybel: a Python wrapper for the OpenBabel cheminformatics toolkit. Chem. Cent. J. 2008, 2, 5.](https://ccj.springeropen.com/articles/10.1186/1752-153X-2-5)
