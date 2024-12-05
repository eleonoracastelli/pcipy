#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:35:47 2024

@author: E Castelli, building on Q Baghi's setup.py file for pylisa
"""

from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open


NAME = "PCIpy"
VERSION = '0.0.1' 
DESCRIPTION = 'Principal Component Interferometry package to process and \
    analyse LISA data.'
# Get the long description from the README file    
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name=NAME, 
        version=VERSION,
        author="Jason Dsouza",
        author_email="<youremail@email.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        keywords='principal component analysis, gravitational-wave interferometry, \
            LISA',
        # python_requires='>=3.5',
        install_requires=['numpy', 'xarray', 'scipy', 'sympy', 
        'pyfftw', 'h5py', 'scikit-learn'],
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'

        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
        ]
)
