# Name: myio/__init__.py
# Purpose: Simple I/O procedure
# Author: Stanislav Stoupin <sstoupin@gmail.com>
#
# Copyright 2012 Argonne National Laboratory
#
# See the file "LICENSE" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.

"""Read and write text files and header.

"""


import os
import re
from numpy import *

#import sys
#import commands
#cmdout=commands.getstatusoutput('echo $HOME')
#libpath0=cmdout[1]+'/bin/'
#sys.path.append(libpath0)
#libpath1=cmdout[1]+'/bin/DTXRD/'
#sys.path.append(libpath1)

#__all__ = ['SffFieldError', 'Sff', 'readSffFile', 'Xafs', 'readXafsFile',
#    'Text', 'readTextFile', 'readFile', 'writeFile']

#class Data(ndarray):
#    """
#    Store columns of data and their associated names.
#
#    These columns can be maniplulated by the get(), ?add()?, and ?set()?
#    methods or by accessing the columns as though they are attributes of the
#    object.
#
#    """
#    def __init__(self, name):
#        """
#        Create an empty dataset.
#        """
#        self.name = name
#        self.__columns = {}
#        self.__columnNames = []
                
#    def get(self, arg):
#        """
#        Retrieve a column's data.        
#        """
#        try: 
#            index=int(arg)-1
#        except ValueError:
#            return None
#        else: 
#            matrix=self.data
#            return matrix[:,index]
        

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
        if char1=='#' or char1=='!' or char1==';' or char1=='{' or char1=='}':
            header=header+line
        elif not(line):
            break
        else:
            line1=line.split()
            if line1!=[] and len(line1)>1: data=data+[line1] 
                        
    data2=array(data, dtype=float)                
    return header, data2

def readFile1(fileName):
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
        if char1=='#' or char1=='!' or char1==';' or char1=='{' or char1=='}' or char1==' ' or char1=='\t':
            header=header+line
        elif not(line):
            break
        else:
            line1=line.split()
            if line1!=[] and len(line1)>1: data=data+[line1] 
                        
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

