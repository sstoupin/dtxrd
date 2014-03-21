#!/usr/bin/env python

from numpy import *

def Fh_dia(h,k,l):
      
      H=2.0*pi*array([h,k,l])
      
      r=[[]]*8
      r[0]=array([0.00,0.00,0.00])
      r[1]=array([0.25,0.25,0.25])
      r[2]=array([0.50,0.50,0.00])
      r[3]=array([0.75,0.75,0.25])
      r[4]=array([0.50,0.00,0.50])
      r[5]=array([0.00,0.50,0.50])
      r[6]=array([0.75,0.25,0.75])
      r[7]=array([0.25,0.75,0.75])
            
      Fh=0
      for x in r:
            Fh=Fh+exp(1j*sum(H*x))
      
      return Fh

print Fh_dia(0,0,0)                              
print Fh_dia(0,0,4)


      
      