#!/usr/bin/env python

'''
a program to calculate x-ray flux based on a scaler counts, detection sensitivity and 
a calibration file of a PIN detector

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys
from numpy import *
#from myio import readFile
from dtxrd.myio import *
import scipy
from scipy.interpolate import interp1d

__version__ = '0.2'

def fatalError(msg):
	sys.stderr.write('Error: ')
	sys.stderr.write(str(msg))
	sys.stderr.write('\n')
	sys.exit(1)

def fatalIOError(err):
	if issubclass(err.__class__, IOError) and err.strerror and err.filename:
		err = '%s: %s' % (err.strerror, err.filename)
	fatalError(err)

def ParseArguments(args):
	try:
		from optik import OptionParser
	except ImportError:
		try:
			from optparse import OptionParser
		except ImportError:
			fatalError(
'This program requires Optik, availible from http://optik.sourceforge.net/\n')

	USAGE = '%prog  COUNTRATE SENSITIVITY_[mA/V] ENERGY_[keV] pin_calibration_file'
	VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>'
	parser = OptionParser(usage=USAGE, version=VERSION)
	opts, args = parser.parse_args(args)

	# did we get the right number of arguments?
	if len(args) != 4:
		parser.print_usage()
		sys.exit(1)

	try:    
	        for k in [0,1,2]:
         	  test = float(args[k])
	except ValueError:
		fatalError('parameters must be valid numbers')

	return opts, args

#wavelength(Angstrom) is calculated and reported according to:
#Energy=h_plank*frequency=momentum*c=h_bar*k*c=2*pi*h_bar*c/lamda

def main(args):

        fn=args[3]
        d1,d2=readFile(fn)
        en=d2[:,0]  #; en_l=list(en)
        cal=d2[:,1] #; cal_l=list(cal)
        
        CR=float(args[0])
        sens=float(args[1])
        en0=float(args[2])
#       
        f_=interp1d(en,cal)
        f=f_(en0)
        print "f = ", f , " [THz/mA]"           
        flux=(CR*1e-5)*f*sens*1.0e12  # [Hz]
        print "flux = ", '%e' %flux, " Hz"
                        
#        if flag==0:
#                print "Energy must be in 5-25 keV interval"                
                
if __name__ == '__main__':
	options, args = ParseArguments(sys.argv[1:])
	main(args)
