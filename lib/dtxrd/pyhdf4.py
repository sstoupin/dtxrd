#!/usr/bin/env python
# program to read hdf4 file
#

import os
from numpy import *

def read_hdf4(fileName):
 
 m=os.popen('hdp dumpsds -n rotation_motor '+fileName)
 stuff0=m.readlines()
 th=float(stuff0[len(stuff0)-2])
      
 f=os.popen('hdp dumpsds -n data '+fileName)
 stuff=f.readlines()

 for x in stuff:
       line=x.rsplit(' ')
       if len(line)>1:
         if line[1]=='Dim0:':
            i_y=stuff.index(x)+1
         elif line[1]=='Dim1:':
            i_x=stuff.index(x)+1   
         elif line[1]=='Data':
#      if x=='\t Data : \n':
            i_d=stuff.index(x)+1
            # print i_d

 linex=stuff[i_x]; linex=linex.split(' '); nx=int(linex[3]); #print 'nx = ', nx
 liney=stuff[i_y]; liney=liney.split(' '); ny=int(liney[3]); #print 'ny = ', ny

 data0=stuff[i_d:len(stuff)]

 data1=[]

 for x in data0:
     x=str(x).split(' ')
     x=x[16:len(x)-1]
     data1=data1+x

 data2=[]
 for k in range(0,ny):
     datax=data1[k*nx:(k+1)*nx]
     data2=data2+[datax]

#########################################################################
# this be better extracted from hdf4 file template
#########################################################################
# dx, dy = 0.06, 0.06  # pixel size [mm]
#########################################################################
 im=array(data2, dtype=float) 
 #print im.dtype
 
 return [th,[nx,ny],im]


