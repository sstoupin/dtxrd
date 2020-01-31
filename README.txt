python-dtxrd - Tools for X-ray diffraction evaluation of single crystals
------------------------------------------------------------------------
dtxrd - calculates reflectivity and transmissivity for a given Bragg
reflection of various single crystals: C (diamond), Si, Ge and Al2O3
(sapphire).

throughput - calculates throughput for a multi-crystal configuration
described by an input file with a special format

rctopo - recombines 2D images of crystal reflectivity at different 
angles on the crystal rocking curve to produce rocking curve images 
(e.g., maps of the rocking curve's width, peak position, etc.)

rctopo-fast - a version of rctopo optimized for fast processing of large datasets (array operations)

seehdf - plots maps of 2D data contained in hdf4 and hdf5 files 

rcpeak - a tool to plot and perform fitting and analysis of a rocking curve
peak data given by a multi-column ASCII file.

specscan - extracts a particular scan from a SPEC file

flux - x-ray flux calculator based on scaler count rate and a detector response

crl - compound reflractive lens calculator

REQUIREMENTS
------------
* hdf (hdf4, includes *hdp* command-line tool)
* python-numpy
* python-scipy
* python-matplotlib
* h5py
* rctopo-fast requires numpy >= 1.17.4

INSTALLATION
------------

The Python Distutils system provides packaging, compilation, and installation
for python-dtxrd

To install, execute the following command as superuser:
  # python setup.py install [OPTIONS]

For more information about installation options, execute the following
command:
  > python setup.py install --help

For information about other Distutils commands, execute the following command:
  > python setup.py --help-commands


AVAILABILITY
------------

* http://python-dtxrd.readthedocs.org
* https://github.com/sstoupin/dtxrd

AUTHOR
------

python-dtxrd was written by Stanislav Stoupin <sstoupin@gmail.com>


Copyright (c) 2014, UChicago Argonne, LLC

All Rights Reserved

XSD, Advanced Photon Source, Argonne National Laboratory


OPEN SOURCE LICENSE

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.  Software changes, 
   modifications, or derivative works, should be noted with comments and 
   the author and organization's name.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation 
   and/or other materials provided with the distribution.

3. Neither the names of UChicago Argonne, LLC or the Department of Energy 
   nor the names of its contributors may be used to endorse or promote 
   products derived from this software without specific prior written 
   permission.

4. The software and the end-user documentation included with the 
   redistribution, if any, must include the following acknowledgment:

   "This product includes software produced by UChicago Argonne, LLC 
   under Contract No. DE-AC02-06CH11357 with the Department of Energy."

****************************************************************************

DISCLAIMER

THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND.

Neither the United States GOVERNMENT, nor the United States Department 
of Energy, NOR uchicago argonne, LLC, nor any of their employees, makes 
any warranty, express or implied, or assumes any legal liability or 
responsibility for the accuracy, completeness, or usefulness of any 
information, data, apparatus, product, or process disclosed, or 
represents that its use would not infringe privately owned rights.

****************************************************************************
