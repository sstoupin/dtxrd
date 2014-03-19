# Name: myio/__init__.py

'''
simple I/O procedure

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

#import os
#import re
from numpy import *

#import sys
#import commands
#cmdout=commands.getstatusoutput('echo $HOME')
#libpath0=cmdout[1]+'/bin/'
#sys.path.append(libpath0)
#libpath1=cmdout[1]+'/bin/DTXRD/'
#sys.path.append(libpath1)
        
def readFile(fileName):
    """Read experiment data and metadata from a file.

    @param fileName: file to open

    @return: header and data
    @rtype:  C{2-tuple} of the form C{(header, data)}

    @raise IOError: could not guess the file format, or open, read, or parse
                    the file
    """    
    header=''
    data=[]    
    input = open(fileName, 'r')
    while 1:
        line=input.readline()
        char1=line[0:1]
        if char1=='#' or char1=='!' or char1==';':
            header=header+line
        elif not(line):
            break
        else:
            line1=line.split()
            if line1!=[]: data=data+[line1] 
                        
    data2=array(data, dtype=float)                
    return header, data2

def readCrystal(fileName):
    header=''
    data=[]    
    input = open(fileName, 'r')
    while 1:
        line=input.readline()
        char1=line[0:1]
        if char1=='#' or char1=='!':
            header=header+line
        elif not(line):
            break        
        else:
            line1=line.split()
            data=data+[line1] 
    return header, data

def coljoin(*args):
    data=[]
    for c in args:
        data=data+[c]
    data2=transpose(array(data))
    return data2

def writeFile(fileName,header,*args):
    dataFile = open(fileName, 'w')
    dataFile.write(header)
    stuff=coljoin(args)
    for x in stuff:
      line=''
      for y in x:
         line=line+' '+"%1.8e"%y
      dataFile.write(line[1:]+'\n')
    dataFile.close

