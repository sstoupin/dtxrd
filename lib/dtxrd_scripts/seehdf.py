#!/usr/bin/env python

'''
a simple hdf4 image viewer

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys
from numpy import *
from pylab import *

import os
if os.path.abspath(os.path.dirname(__file__)).split(os.sep)[-2] == 'lib':
    '''when running this script from the source directory'''
    sys.path.insert(0, os.path.abspath('..'))

from dtxrd.pyhdf4 import *
from matplotlib.patches import Patch
############################################################
dx, dy = 0.06, 0.06  # CCD camera pixel size [mm]
############################################################
__version__='0.11'
proginfo = 'viewhdf' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>'


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

        USAGE = '%prog [OPTIONS...] filename1 filename2 ...'
        VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>\n' 
        parser = OptionParser(usage=USAGE, version=VERSION)                        
        #
        # Behavior
        #
        parser.add_option('-o', '--output',
                action='store',
                dest='output',
                default=None,
                help='write results to file F (defaults to stdout)',
                metavar='F')

        opts, args = parser.parse_args(args)
        if len(args) < 1:
                parser.print_usage()
                sys.exit(1)
        return opts, args  


def main():
  opts, args = ParseArguments(sys.argv[1:])
  if opts.output is not None:
     try:
        outFile = open(opts.output, 'w')
     except IOError, e:
        fatalIOError(e)
  else:
     outFile = sys.stdout
  
  count=1
  for fn in args:   
    angle, size, im = read_hdf4(fn)
    nx=size[0] ;print nx
    ny=size[1] ;print ny
    
    im=fliplr(im)
    xyrange=(0,dx*nx,0,dy*ny)
    plt.figure(count)
    plt.imshow(im, extent=xyrange)  
    plt.colorbar()
    plt.title(fn)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]') 
    count=count+1  
  
  outFile.write('# seehdf '+__version__+' by Stanislav Stoupin <sstoupin@gmail.com>\n')
  outFile.write('# image: '+str(nx)+' x '+str(ny)+'\n')
  outFile.close  
  plt.show()  


if __name__ == '__main__':
    main()
