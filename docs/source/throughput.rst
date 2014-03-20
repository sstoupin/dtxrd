
.. _throughput:

************
throughput
************

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

calculate throughput of a multicrystal configuration 
(dynamical theory of x-ray diffraction for perfect crystals)

SYNOPSIS
============

       throughput [options] func dpsi Ec dEx ne input_file

DESCRIPTION
============

A program to calculate throughput of a multicrystal configuration given by input_file using dynamical theory of
x-ray diffraction for perfect crystals
For a brief summary of options and parameters run::

    throughput -h

PARAMETERS
============

:func: angular divergence distribution function for incoming x-rays:
       
       ====  ============  ==================
       func  description   dpsi meaning
       ====  ============  ==================
       g     Gaussian      rms
       l     Lorentzian    fwhm
       ====  ============  ==================

:dpsi: angular divergence in units of [urad]
       
       * rms for Gaussian distribution
       * fwhm for Lorentzian distribution

:Ec:   central  energy  for  energy  distribution  in  units of [keV] (corrected automatically depending on the
       choice of source and presence of crystals in backscattering configuration in the input file)

:dEx:    energy half-range in units of [meV]

:ne:     number of steps in the energy grid to perform the calculation (max 1000)

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

:-a NTH, --angular_scan=NTH:
       perform crystal rotation (angular scan) with NTH points

:-s SRC, --source=SRC:
       type of energy distribution for the source:
       
       ====    ===============================================
       SRC     description
       ====    ===============================================
       0       flat distribution (default)
       1       Cu K-alpha source
       9       energy distribution from file source_e.dat
       ====    ===============================================

EXAMPLES
===========

This is an example input file for calculation of throughput of a multicrystal configuration

:download:`thru_hhlmC.in <../../examples/throughput/thru_hhlmC.in>`

.. literalinclude:: ../../examples/throughput/thru_hhlmC.in
   :language: guess
   :linenos:

This is an example input file for calculation of a rocking curve in a multicrystal configuration

:download:`c2rc_hhlmC.in <../../examples/throughput/c2rc_hhlmC.in>`

.. literalinclude:: ../../examples/throughput/c2rc_hhlmC.in
   :language: guess
   :linenos:


SEE ALSO
============

* :ref:`dtxrd`
* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
