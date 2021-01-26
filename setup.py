#External package imports 

from setuptools import setup

#Creating a valid setup.py file

setup(name = 'mbfit',
	  description = '''A package that enables you to create, and visualize 1 body and 2 body MB
	                   and TTM fits.''', 
	  author = 'Paesani Lab, University of California, San Diego',  # we should add a comma seperated list of authors,
                                                                    # but I am not sure of the entire list and I don't
                                                                    # want to leave anyone out so I'll wait to do so
                                                                    # until I know who to include/
	  license = 'University of California San Diego',
	  author_email = 'xyz@abc.com', #Will update this later!
	  packages=['mbfit', 'mbfit.utils', 
	            'mbfit.fitting', 'mbfit.configurations',
	            'mbfit.database', 'mbfit.polynomials',  
	            'mbfit.exceptions', 'mbfit.molecule', 
	            'mbfit.calculator' ],
	  zip_safe = False
	)
