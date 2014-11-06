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
   'numpy>=1.4.1', 
   'scipy>=0.7.2', 
   'matplotlib>=0.99.1.1',
]

setup(
    name = 'python-dtxrd',
    version = '1.1',
    url = 'https://subversion.xray.aps.anl.gov/dtxrd/',
    maintainer = 'Stanislav Stoupin',
    maintainer_email = 'sstoupin@gmail.com',
    license = 'UChicago Argonne, LLC OPEN SOURCE LICENSE',
    description = 'Tools for X-ray diffraction evaluation of single crystals',
    install_requires = install_requires,
    package_dir = {'': 'lib'},
    #packages = ['dtxrd', 'dtxrd.myio'],
    packages = find_packages('lib'),
    package_data = {'dtxrd': ['asf/*.asf', 'data/*.dat']},
    scripts = ['dtxrd', 'throughput', 'rctopo', 'seehdf', 'rcpeak', 'specscan', 'flux'],
    zip_safe = False,
)
