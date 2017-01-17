
.. _rctopo:

************
rctopo
************

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

x-ray rocking curve topography program 

SYNOPSIS
============

::

       rctopo [options] filename1 filename2 ... filenameN

DESCRIPTION
============

A program to process a sequence of topographs collected at different angles on the 
rocking curve of a crystal to generate maps of the rocking curve parameters.
Supported area detector file formats: HDF4 (.hdf) and HDF5 (.h5)

OPTIONS
============

For a brief summary run::

    rctopo -h

:-v,   --version:
       show program's version

:-h,         --help:
       show summary of options

:-o F, --output=F:
       write calculated results to file F (default to stdout); also, generates output 1pix_total.dat
       showing the rocking curve from the central pixel and the total rocking curve

:-w D, --output=D:
       write slice data to file D (default: do not write)

:-t T, --threshold=T:
       threshold for data processing to reject "weak" rocking curves to define
       crystal boundaries (default T=1.05)

:-b bkg, --background=bkg:
       user defined background (dark current) of the area detector (default value is estimated
       from the rocking curve tails)

:-r STRING, --range=STRING:
       xy-range for display and analysis (STRING='x1 x2 y1 y2', where x1,x2,y1,y2 are in units of
       [mm])

:-x CONST, --xslice=CONST:
       slice and plot distributions at a fixed coordinate X = CONST

:-y CONST, --yslice=CONST:
       slice and plot distributions at a fixed coordinate Y = CONST

:-f CONST, --factor=CONST:
       scale colormap range on topographs by CONST*FWHM_av, where FWHM_av is the average FWHM

:-m CONST, --magnify=CONST:
       shrink range on the colormaps by factor CONST; applies only to the colormaps which 
       represent characteristic width of the rocking curve width (FWHM and STDEV)

:-n STRING, --name=STRING:
       include sample name STRING in the figure title

:-d CONST, --deglitch=CONST:
       deglitch rocking curve data with CONST as a threshold parameter (e.g., CONST=1.1) (default - no deglitching)

:-g,   --gaussian:
       perform Gaussian curve fitting (smoothes noisy images)

:-s,   --transpose:
       transpose image array for plotting

:-u uname, --units=uname:
       assign the original angular units (uname): deg, arcsec or urad (default - deg)

:-p,   --publish:
       generate additional figures with publication quality (requires a separate figures.py script)

:-c,   --conduct:
       process forward diffraction data       

:-i,   --instrument:      
       read the detector and the image analysis parameters from an instrument file ccd.py

:-z CONST, --integrate=CONST:
       integrate reflectivity and normalize by the theoretical angular acceptance for perfect crystal (CONST);
       specific cases: z = 0 - no integration (default), z = -1 - integrate and normalize by the maximum value 

:-e SPECSCAN, --external=SPECSCAN:
	read angular steps from the first column of a SPEC scan 

EXAMPLES
===========

This archive contains a set of hdf images of a diamond 111 crystal collected at 
different angles on its rocking curve using a Cu :math:`K_{\alpha}` radiation collimated by a 
strongly asymmetric Si 220 reflection. 

:download:`SA1.zip <../../examples/rctopo/SA1.zip>`

to perform quick evaluation run::

    rctopo -s -u deg *hdf

.. image:: ../../examples/snapshots/rctopo00.png
            :width: 50 %
	    :alt: diamond SA1 	    	    

to better define crystal boundary (threshold for analysis) and to obtain a smooth image (Gaussian fitting for each pixel) run::

    rctopo -t 1.1 -g -s -u deg *hdf

.. image:: ../../examples/snapshots/rctopo0.png
            :width: 50 %
	    :alt: diamond SA1 fitting/threshold

to display the name of the sample in the figure title run::

    rctopo -t 1.1 -g -s -u deg -n diamond1 *hdf

.. image:: ../../examples/snapshots/rctopo1.png
            :width: 50 %
	    :alt: diamond SA1 name	    	    

to perform statistical analysis and visualization over a specified region run::

    rctopo -r '1.5 3.5 4 6' -t 1.1 -g -s -u deg -n diamond1 *hdf

.. image:: ../../examples/snapshots/rctopo2.png
            :width: 50 %
	    :alt: diamond SA1 working region	    	    


SEE ALSO
============

* :ref:`seehdf`
* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
