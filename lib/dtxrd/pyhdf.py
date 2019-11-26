#!/usr/bin/env python
# a subroutine to read hdf4 and hdf5 files
#
# 11/12/2019 
# replaced deprecated ".value" method to [()] for h5py 

import os
import h5py
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
    im = flipud(im)
 
    return [th,[nx,ny],im]

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def read_hdf5(fileName,rbin,data_path,th_path,chi_path):

    hdf5 = h5py.File(fileName, 'r') 

    # EXTRACT DATA    
    try:
        ds = hdf5[data_path]
    except:
        msg = 'datapath not found in HDF5 file:'
        msg += '\n  file: ' + fileName
        msg += '\n  path: ' + data_path
        raise IOError(msg)
    #image_data = numpy.ma.masked_less_equal(ds.value, bkg)
    #image_data = image_data.filled(image_data.min())                                              
    #im = ds.value
    im = ds[()]
    im = array(im, dtype=float)    
    shape0 = ds.shape
    if len(shape0) == 3:
        im=im[0]    
        (n0,ny,nx) = shape0
    elif len(shape0) == 2:
        (ny,nx) = shape0
    else: 
        msg = 'data of unknown structure:'
        msg +='\n file: ' +fileName        
        msg += '\n  path: ' + data_path
        raise IOError(msg)

    # REBIN DATA                         
    #print("nx = ", nx)
    #print("ny = ", ny)
    # rebin 
    shape1 = (ny//rbin,nx//rbin)
    im = rebin(im,shape1)
    (ny,nx) = shape1
    #im = flipud(im)
                
    # EXTRACT ANGLES
    try:
        thentry = hdf5[th_path]
    except:
        msg = 'datapath not found in HDF5 file:'
        msg += '\n  file: ' + fileName
        msg += '\n  path: ' + th_path
        raise IOError(msg)
    #th = float(thentry.value)
    th = float(thentry[()])
    #
    try:
        chientry = hdf5[chi_path]
    except:
        msg = 'datapath not found in HDF5 file:'
        msg += '\n  file: ' + fileName
        msg += '\n  path: ' + chi_path
        raise IOError(msg)
    #chi = float(chientry.value)
    chi = float(chientry[()])
    hdf5.close()                     
    return [th,chi,[nx,ny],im]
    
