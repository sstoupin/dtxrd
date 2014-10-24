#!/usr/bin/env python

'''
a program to calculate x-ray flux based on a scaler counts, detection sensitivity and 
a calibration file of a PIN detector

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''
import os
import sys
from numpy import *
from myio import readFile
import scipy
from scipy.interpolate import interp1d

prog = os.path.basename(sys.argv[0])        
__version__ = '0.21'

def fatalError(msg):
	sys.stderr.write('Error: ')
	sys.stderr.write(str(msg))
	sys.stderr.write('\n')
	sys.exit(1)

def fatalIOError(err):
	if issubclass(err.__class__, IOError) and err.strerror and err.filename:
		err = '%s: %s' % (err.strerror, err.filename)
	fatalError(err)

def ParseArguments():
        import argparse   # requires Python 2.7 or higher
        msg = prog + '  version: ' + __version__ + '\n'*2 + __doc__.strip()+'\n'
        msg1 = 'countrate from a scaler [counts/s] \n'
        msg2 = 'pre-amp sensitivity [mA/V] \n'
        msg3 = 'photon energy [keV] \n'
        msg4 = 'pin detector calibration file containing two columns: 1 - energy[keV] 2 - response[THz/mA] \n'
        msg5 = 'write calculated parameters to file F (defaults to stdout) \n'
        #
        parser = argparse.ArgumentParser(prog=prog, description=msg, formatter_class = argparse.RawDescriptionHelpFormatter) 
        parser.add_argument('-v', '--version', action='version', version=__version__)
        #                
        parser.add_argument('CR', action='store', type=float, nargs=1, help=msg1)
        parser.add_argument('sens', action='store', type=float, nargs=1, help=msg2)
        parser.add_argument('en0', action='store', type=float, nargs=1, help=msg3)
        parser.add_argument('fn', action='store', help=msg4)
        parser.add_argument('-o', '--output', action='store', dest='output', default=None, help=msg5, metavar='F')
        #        
        return parser.parse_args()
        
def main():
        stuff = ParseArguments()
        
        if stuff.output is not None:
            try:
               outFile = open(cmd_opts.output, 'w')
            except IOError, e:
               fatalIOError(e)
        else:
            outFile = sys.stdout
########################################################################################                    
        d1,d2=readFile(stuff.fn)
        en=d2[:,0]  #; en_l=list(en)
        cal=d2[:,1] #; cal_l=list(cal)
                
        CR=float(stuff.CR[0])
        sens=float(stuff.sens[0])
        en0=float(stuff.en0[0])
#       
        r_=interp1d(en,cal)
        R=r_(en0)
        flux=(CR*1e-5)*R*sens*1.0e12  # [Hz]
        ################################################################################
        ## OUTPUT      
        ################################################################################  
        outFile.write('##############################################################\n')
        outFile.write('##### ' + prog + ' v'+__version__+' ############################################\n')
        outFile.write('##### Author: Stanislav Stoupin ## sstoupin@aps.anl.gov ######\n')
        outFile.write('##############################################################\n')                                
        outFile.write('Ex [keV] = ' + str(en0) +' photon energy \n')
        outFile.write('R [THz/mA] = '+ str(R) + ' pin detector response \n')        
        outFile.write('flux [Hz] = ' + '%e' %flux + ' photon flux \n')        
                
if __name__ == '__main__':
	main()
