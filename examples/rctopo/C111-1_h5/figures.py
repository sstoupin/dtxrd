from pylab import *
from constants import r2d
r2d = 180.0/pi
import numpy as np
def figplot(s,dx,dy,indx1,indx2,indy1,indy2,peak,fwhm,stdev,com,thmid,thneg,thpos):
##
  thc=21.96
  thsi=23.6449
  dth = 2.0*(thsi-thc)
  dy=dy/(sin(thc/r2d+dth/r2d) - cos(thc/r2d+dth/r2d)*tan(dth/r2d))
  #dx1=dx
  #if s==1:
  xyrange=(0.0,dx*(indx2-indx1),0.0,dy*(indy2-indy1))
  #else:
  #   xyrange=(0.0,dy*(indy2-indy1),0.0,dx1*(indx2-indx1))
                                  
  matplotlib.rcParams.update({'font.size': 20})
  matplotlib.rcParams['xtick.major.pad']='8'
  matplotlib.rcParams['ytick.major.pad']='8'
  
  f31=plt.figure(9)
  shr1=0.8
  pad1=0.02
  asp1=20.0
  lbp1=20.0
  xti = [0,1,2,3,4,5,6,7,8,9,10,11]
  #xti = [0,1,2,3]
  yti = [0,1,2,3,4,5,6,7,8,9,10,11]
  #yti = [0,1,2,3,4,5,6,7]
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
  #--------------------------------------------------------------------------------------
  plt.subplot(141)
  vmin0=26.5
  vmax0=29.5
  trmax=10.0*vmax0
  fwhm_m=np.ma.array(fwhm,mask=peak>trmax)
  imgplot = plt.imshow(fwhm_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
# imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[26.5,27.0,27.5,28.0,28.5,29.0,29.5],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('FWHM')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(142)
  vmin0=-2.0
  vmax0=2.0
  trmax=10.0*vmax0
  thmid_m=np.ma.array(thmid,mask=peak>trmax)
  imgplot = plt.imshow(thmid_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('mid-point')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(143)
  vmin0=-2.0
  vmax0=2.0
  trmax=10.0*vmax0
  thneg_m=np.ma.array(thneg,mask=peak>trmax)
  imgplot = plt.imshow(thneg_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('left slope')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(144)
  vmin0=-2.0
  vmax0=2.0
  trmax=10.0*vmax0
  thpos_m=np.ma.array(thpos,mask=peak>trmax)
  imgplot = plt.imshow(thpos_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('right slope')
  plt.xticks(xti)
  plt.yticks(yti)




     