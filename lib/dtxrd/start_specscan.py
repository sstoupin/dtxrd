#!/usr/bin/env python

'''
extracts scans from a SPEC file returning individual files with headers

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import os
import sys
#from numpy import *

__version__='0.02'
proginfo = 'specscan v' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>'

########################################################################
# modifications:
########################################################################
#v0.02  now writes data only upon request (-w option)
#       by default only prints columns
#---------------------------------------------------
# Functions
#---------------------------------------------------

cbank=('b','g','r','c','m','y')

#############################################################################################

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
    import argparse
    prog = os.path.basename(sys.argv[0])
    msg = prog + '  version: ' + __version__ + '\n'*2 + __doc__.strip()
    parser = argparse.ArgumentParser(prog=prog, 
                                     description=msg,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', 
                        action='version', 
                        version=__version__)

    # Output
    msg = 'write extracted scans into files <scan_number>.dat in current directory'
    msg += ' (default: print scan command and column names but do not write output files)'
    parser.add_argument('-w', '--write',
            action='store_true',
            dest='write',
            default=False,
            help=msg)

    # Positional Parameters
    parser.add_argument('scan_number', 
                        action='store', 
                        type=int,
                        nargs='+', 
                        help="one or more scan numbers")
    parser.add_argument('fileName', 
                        action='store', 
                        help="SPEC data file name")
    return parser.parse_args()


##################################################################################################
## Main stuff
##################################################################################################


def main():
#     options, args = ParseArguments(sys.argv[1:])
    cmd_opts = ParseArguments()

    #--------------------------------------------------------
    # initialize variables :
    #--------------------------------------------------------

    count=0
    #--------------------------------------------------------
    # do stuff
    #--------------------------------------------------------

    fileName = cmd_opts.fileName
    scans    = cmd_opts.scan_number

    for fn in scans:
        raw=open(fileName, 'r')
        count=count+1
        #
        header=''
        data=''
        while 1:
            line=raw.readline()
            if not(line):
                break
            if line[0:2]=='#S':
                splitline=line.split()
                if splitline[1]==str(fn):
                    print line
                    header=header+line
                    #
                    while 2:
                        line1=raw.readline()
                        char1=line1[0:1]
                        lst1=line1.split()
                        if lst1==[]:
                            break
                        elif not(line1):
                            break
                        elif char1=='#':
                            header=header+line1
                            if lst1[0]=='#L':
                                print 'columns in the scan:'
                                names=lst1[1:]
                                numbers=range(1,len(names)+1)
                                for (i,x) in map(None,numbers,names):
                                    print 'col '+str(i)+': '+x
                        else:
                            data=data+line1

        if cmd_opts.write:
            fname = str(fn)+'.dat'
            outFile=open(fname, 'w'); print "writing extracted scan into file: " + fname
            outFile.write(header)
            outFile.write(data)
            outFile.close()
            raw.close()


if __name__ == '__main__':
    main()
