#!/usr/bin/env python

# Name: setup.py
# Purpose: python-dtxrd distutils install program
# Author: Stanislav Stoupin <sstoupin@gmail.com>
#
# Copyright 2014
#
# See the file "LICENSE" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.


from setuptools import setup, find_packages

install_requires = [
   'numpy>=1.12.1', 
   'scipy>=0.18.1',
   'matplotlib>=0.99.1.1',
   'h5py>=2.0.1',
   'pillow>=4.0.0',
]

setup(
    name = 'python-dtxrd',
    version = '1.8',
    url = 'https://github.com/sstoupin/dtxrd',
    maintainer = 'Stanislav Stoupin',
    maintainer_email = 'sstoupin@gmail.com',
    license = 'UChicago Argonne, LLC OPEN SOURCE LICENSE',
    description = 'Tools for X-ray diffraction evaluation of single crystals',
    install_requires = install_requires,
    package_dir = {'': 'lib'},
    #packages = ['dtxrd', 'dtxrd.myio'],
    packages = find_packages('lib'),
    package_data = {'dtxrd': ['asf/*.asf', 'data/*.dat']},
    scripts = ['dtxrd', 'throughput', 'rctopo', 'rctopo-fast', 'seehdf', 'rcpeak', 'specscan', 'flux', 'crl'],
    zip_safe = False,
)
