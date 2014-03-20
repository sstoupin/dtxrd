
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

:phi:    azimuthal angle of incidence

:T:      crystal temperature [K]

:d:      crystal thickness [mm]

:flag: =====   =================================================================
       flag    description
       =====   =================================================================
       a       perform calculation at a given glancing angle of incidence theta
       e       perform calculation at a given photon energy Ex
       =====   =================================================================

:theta: glancing angle of incidence, theta (:math:`\theta`)

:Ex: photon energy, Ex (:math:`E_{\mathrm X}`)


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

EXAMPLES
===========

to calculate a rocking curve of a 1-mm-thick Si (111) crystal at 8 keV (111 reflection, Bragg case) run::

       dtxrd Si 1 1 1 0 0 300 1 e 8

.. image:: ../../examples/snapshots/Si111_8keV.png
            :width: 90 %
	    :alt: Si111 at 8keV

to calculate a rocking curve of a 0.1-mm-thick C (001) crystal at 12 keV (220 reflection, Laue case) run::

       dtxrd C 2 2 0 45 0 300 0.1 e 12 

.. image:: ../../examples/snapshots/C220_Laue.png
            :width: 90 %
	    :alt: C220 Laue at 12keV


SEE ALSO
============

* :ref:`throughput`
* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
