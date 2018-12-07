#External package imports 

from setuptools import setup

#Creating a valid setup.py file

setup(name = 'potential_fitting',
	  description = '''A package that enables you to create, and visualize 1 body and 2 body MB
	                   and TTM fits.''', 
	  author = 'Paesani Lab',
	  license = 'University of California San Diego',
	  author_email = 'xyz@abc.com', #Will update this later!
	  packages=['potential_fitting', 'potential_fitting.utils', 
	            'potential_fitting.fitting', 'potential_fitting.configurations',
	            'potential_fitting.database', 'potential_fitting.polynomials',  
	            'potential_fitting.exceptions', 'potential_fitting.molecule', 
	            'potential_fitting.calculator' ],
	  zip_safe = False
	)
