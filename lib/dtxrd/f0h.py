#!/usr/bin/env python

from numpy import *
from scipy import interpolate


hpl=4.13566733e-15 # [eV*s] Planck constant
c=299792458e+10    # [Angstrom/s] speed of light

def f0h_ICD(atom,th,Ex):
    lamx=hpl*c/Ex   #; print "lamx [A] = ", lamx
    qx=sin(th)/lamx #; print "qx [A-1] = ", qx
        
    fileName='f0h_ICD.dat'
    data=open(fileName, 'r')
    while 1:
        line=data.readline()
        stuff=line.split(' ')
        if stuff[0]=='Q':
          line1=data.readline()
          q=line1.split(' ')
        elif stuff[0]==atom:
          line1=data.readline()
          f0=line1.split(' ')
          break
        elif not(line): 
          break
    
    q=array(q, dtype=float)   #; print q
    f0=array(f0, dtype=float) #; print f0
    f0q=interpolate.interp1d(q,f0)
    return f0q(qx)

#from f_cromer import *
#print f_cromer('Si',10080)
#print f0h('Si',0.0,10080)

