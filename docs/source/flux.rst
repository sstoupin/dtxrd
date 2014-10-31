
.. _flux:

************
flux
************

x-ray flux calculator

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

SYNOPSIS
============

::

       flux [options] CR sens Ex pinXXX.dat

DESCRIPTION
============

A program to calculate x-ray flux based on a scaler counts, detection sensitivity and 
a calibration file of a PIN detector. A gain 1X is assumed on the voltage-to-frequency converter.


INPUT PARAMETERS
=================

:CR: count rate from a scaler [counts/s]

:sens: sensitivity of a current amplifier [mA/V]

:Ex: photon energy [keV]

:pinXXX.dat: detector calibration file containing two columns: 

             ==================  ===========================
             photon energy[keV]  detector response [THz/mA]
             ==================  ===========================
    
OPTIONS
============

:-v, --version:
       show program's version number

:-h, --help:
       show summary of options.

:-o F, --output=F:
       write results to file F (default to stdout)

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
