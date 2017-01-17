
.. _crl:

************
crl
************

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

x-ray refractive lens gain calculator

SYNOPSIS
============

::

       crl [options] element S1 S2 sh sv R0 da N Ex sig

DESCRIPTION
============

A program to calculate parameters of a compound refractive lens 


INPUT PARAMETERS
=================

:element:
       material: C (diamond), Si (silicon), Ge (germanium) or Al2O3 (sapphire), etc.

:S1:
       distance from the source to the lens [m]

:S2:
       distance from the lens to the imaging plane [m]    

:sh: 
       horizontal source size FWHM [um]

:sv:
       vertical source size FWHM [um]

:R0:
       half-aperture of the lens [mm] 

:da:  
       distance between apexes of the lens [um]

:N:    
       number of lenses in the stack [integer] 

:Ex:   
       photon energy [keV]sig - rms surface roughness [um]

:sig:
       r.m.s. surface roughness [um]


OPTIONS
============

:-v, --version:
       show program's version number and exit

:-h, --help:
       show the help message

:-o F, --output=F:
       write results of calculations to file F (default to stdout)

:-w D, --write=D:
       write profile data to file D (default - no action)

SEE ALSO
============

* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
