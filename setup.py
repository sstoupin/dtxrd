#!/usr/bin/env python

# Name: setup.py
# Purpose: python-dtxrd distutils install program
# Author: Stanislav Stoupin <sstoupin@gmail.com>
#
# Copyright 2014
#
# See the file "LICENSE" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.


import distutils
from distutils.core import setup

setup(
    name = 'python-dtxrd',
    version = '1.0',
    url = 'https://subversion.xray.aps.anl.gov/dtxrd/',
    maintainer = 'Stanislav Stoupin',
    maintainer_email = 'sstoupin@gmail.com',
    license = 'UChicago Argonne, LLC OPEN SOURCE LICENSE',
    description = 'Tools for X-ray diffraction evaluation of single crystals',
        package_dir = {'': 'lib'},
        packages = ['dtxrd', 'dtxrd.myio'],
	package_data = {'dtxrd': ['asf/*.asf', 'data/*.dat']},
	scripts = ['dtxrd', 'throughput', 'rctopo', 'seehdf', 'rcpeak', 'specscan']	
)
