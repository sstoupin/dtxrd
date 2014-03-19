
.. _dtxrd:

************
dtxrd
************

x-ray diffraction calculator 
(dynamical theory of x-ray diffraction for perfect crystals)

SYNOPSIS
============

::

       dtxrd [options] crystal h k l eta phi T d ["a" | "e"] [theta | Ex]

For example, a Si (111) rocking curve at 8 keV::

       dtxrd Si 1 1 1 0 0 300 1 e 8

DESCRIPTION
============

A  program to calculate parameters of a single Bragg reflection for 
a monochromatic incident wave using dynamical theory of x-ray diffraction for perfect crystals
For a brief summary run::

    dtxrd -h

PARAMETERS
============

:crystal:
       crystal type: C (diamond), Si (silicon), Ge (germanium) or Al2O3 (sapphire)

:h k l:  Miller indicies of a chosen Bragg reflection

:eta:    asymmetry angle

:phi:    asimuthal angle of incidence

:T:      crystal temperature [K]

:d:      crystal thickness [mm]

:flag: =====   =================================================================
       flag    description
       =====   =================================================================
       a       perform calculation at a given glancing angle of incidence theta
       e       perform calculation at a given photon energy Ex
       =====   =================================================================

:theta: glancing angle of incidence, theta (:math:`\theta`)

:Ex: photon energy, Ex (:math:`E_x`)


OPTIONS
============

:-v, --version:
       show version of program.

:-h, --help:
       show summary of options.

:-o F, --output=F:
       write results to file F (default to stdout)

:-w D, --write=D:
       write data to file D (default - no action)

:-p, --pi:
       :math:`\pi` polarization for incident wave (default - :math:`\sigma` polarization)

:-c, --conv:
       convolve data with a virtual instrumental resolution function having FWHM of 1/10 of  the  Darwin  width
       and report the resulting FWHM of the reflectivity curve

SEE ALSO
============

* :ref:`throughput`
* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
