
.. _rcpeak:

************
rcpeak
************

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

plot and calculate parameters of a reflectivity curve

SYNOPSIS
============

::

       rcpeak [options] filename1 filename2 ... filenameN

DESCRIPTION
============

A  program  to  plot and calculate parameters of a reflectivity curve or any kind of distribution with a single
peak

OPTIONS
============

For a brief summary run::

  rcpeak -h


:-v, --version:
       show version of program.

:-h, --help:
       show summary of options.

:-o F, --output=F:
       write results to file F (default to stdout)

:-p, --plot:
       plot fitted model distributions (Gaussian and Lorentzian)

:-c, --contrast:
       calculate spectral contrast at 1.578*FWHM

:-u unit, --unit=unit:
       conversion factor for x-axis

:-m mult, --mult=mult:
       multiplication factor for y-axis

:-x Nx, --name_x=Nx:
       x-axis label/unit

:-y Ny, --name_y=Ny:
       y-axis label/unit

:-b bkg, --bkg=bkg:
       user defined background for calculation of FWHM
       (by default determined from the peak tails)

:-t n_th, --th=n_th:
       x-values (e.g., rocking curve angle th) are in column number n_th (default n_th=1)

:-r n_r, --r=n_r:
       y-values (e.g., reflectivity values r) are in column number n_r (default n_r=2)

:-d nder, --deriv=nder:
       perform interpolation with number of points, calculate and plot
       numerical derivative of data (e.g., for beam size estimation in a knife-edge scan)

:-n, --normalize:
       normalize data for plotting (e.g., userful for visual comparision of FWHM
       of peaks with different maximum intensity)

SEE ALSO
============

* :ref:`dtxrd`
* :ref:`throughput`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
