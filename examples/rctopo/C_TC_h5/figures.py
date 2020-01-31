from pylab import *
r2d = 180.0/pi
import numpy as np
def figplot(s,dx,dy,indx1,indx2,indy1,indy2,peak,fwhm,stdev,com,thmid,thneg,thpos):
##
  #thc=21.96
  #thsi=23.6449
  #dth = 2.0*(thsi-thc)
  #dy=dy/(sin(thc/r2d+dth/r2d) - cos(thc/r2d+dth/r2d)*tan(dth/r2d))
  #dx1=dx
  dy=dy/cos(18.0/r2d)
  #if s==1:
  xyrange=(0.0,dx*(indx2-indx1),0.0,dy*(indy2-indy1))
  #else:
  #   xyrange=(0.0,dy*(indy2-indy1),0.0,dx1*(indx2-indx1))
                                  
  matplotlib.rcParams.update({'font.size': 16})
  matplotlib.rcParams['xtick.major.pad']='8'
  matplotlib.rcParams['ytick.major.pad']='8'
  
  f31=plt.figure(9)
  shr1=0.8
  pad1=0.02
  asp1=20.0
  lbp1=20.0
  #xti = [0,1,2,3,4,5,6,7,8,9,10,11]
  xti = [0.0,0.2,0.4,0.6,0.8,1.0]
  #yti = [0,1,2,3,4,5,6,7,8,9,10,11]
  yti = [0.00,0.05,0.10,0.15,0.20,0.25]
  #-------------------------------------------------------------------------------------
  #plt.subplot(131)    
  #vmin0=0
  #vmax0=1.0
  #trmax=2.0*vmax0
  #peak_m=np.ma.array(peak,mask=peak>trmax)
  #imgplot = plt.imshow(peak_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  #divider = make_axes_locatable(ax)
  #cax = divider.append_axes("right", size="5%", pad=0.05)  
  ysize=xyrange[3]; print ysize
  xsize=xyrange[1]; print xsize
  print "y/x = ", ysize/xsize
  #plt.colorbar(imgplot,ticks=[0,0.2,0.4,0.6,0.8,1.0],shrink=shr1,pad=pad1,aspect=asp1)
#  plt.colorbar(imgplot,ticks=[-20,-10,0,10,20],orientation='horizontal')
  #plt.xlabel('x [mm]',labelpad=lbp1)
  #plt.ylabel('y [mm]',labelpad=lbp1)
  #plt.title('Intensity [n.u.]')
  #plt.xticks(xti)
  #plt.yticks(yti)
  #----------------------------------------------------------------------------------------------------
  plt.subplot(131)
  vmin0=0
  vmax0=0.3
  trmax=10.0*vmax0
  peak_m=np.ma.array(peak,mask=peak>trmax)
  imgplot = plt.imshow(peak_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
# imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[0.0,0.1,0.2,0.3],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('peak intensity [n.u.]')
  plt.xticks(xti)
  plt.yticks(yti)
  #----------------------------------------------------------------------------------------------------
  plt.subplot(132)
  vmin0=1.5
  vmax0=3.0
  trmax=10.0*vmax0
  fwhm_m=np.ma.array(fwhm,mask=peak>trmax)
  imgplot = plt.imshow(fwhm_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
# imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[1.5,2.0,2.5,3.0],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('peak width [$\mu$rad]')
  plt.xticks(xti)
  plt.yticks(yti)
  #---------------------------------------------------------------------------------------------------- 
  plt.subplot(133)
  vmin0=-1.0
  vmax0=1.0
  trmax=10.0*vmax0
  com_m=np.ma.array(com,mask=peak>trmax)
  imgplot = plt.imshow(com_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-1.0,0.0,1.0],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('peak position (COM) [$\mu$rad]')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 




     
