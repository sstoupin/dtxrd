from pylab import *
import numpy as np
r2d = 180.0/np.pi
def figplot(s,dx,dy,indx1,indx2,indy1,indy2,peak,fwhm,stdev,com,thmid,thneg,thpos):
##
  thc=21.96
  thsi=23.6449
  dth = 2.0*(thsi-thc)
  #dy=dy/(sin(thc/r2d+dth/r2d) - cos(thc/r2d+dth/r2d)*tan(dth/r2d))
  dy = dy*cos(dth/r2d)/sin(thc/r2d)
  #dx1=dx
  #if s==1:
  xyrange=(0.0,dx*(indx2-indx1),0.0,dy*(indy2-indy1))
  #else:
  #   xyrange=(0.0,dy*(indy2-indy1),0.0,dx1*(indx2-indx1))
                                  
  matplotlib.rcParams.update({'font.size': 20})
  matplotlib.rcParams['xtick.major.pad']='8'
  matplotlib.rcParams['ytick.major.pad']='8'
  
  f31=plt.figure(figsize=(16,12))
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
  ysize=xyrange[3] #; print(ysize)
  xsize=xyrange[1] #; print(xsize)
  #print("y/x = ", ysize/xsize)
  #plt.colorbar(imgplot,ticks=[0,0.2,0.4,0.6,0.8,1.0],shrink=shr1,pad=pad1,aspect=asp1)
#  plt.colorbar(imgplot,ticks=[-20,-10,0,10,20],orientation='horizontal')
  #plt.xlabel('x [mm]',labelpad=lbp1)
  #plt.ylabel('y [mm]',labelpad=lbp1)
  #plt.title('Intensity [n.u.]')
  #plt.xticks(xti)
  #plt.yticks(yti)
  #--------------------------------------------------------------------------------------
  plt.subplot(221)
  vmin0=27.0
  vmax0=30.0
  trmin=5500
  trmax=65536
  fwhm_m1=np.ma.array(fwhm,mask=peak<trmin)
  fwhm_m2=np.ma.array(fwhm_m1,mask=peak>trmax)  
  imgplot = plt.imshow(fwhm_m2, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
# imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[27,28,29,30],shrink=shr1,pad=pad1,aspect=asp1)
  #plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('FWHM')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(222)
  vmin0=-2
  vmax0=2
  thmid_m1=np.ma.array(thmid,mask=peak<trmin)
  thmid_m2=np.ma.array(thmid_m1,mask=peak>trmax)  
  imgplot = plt.imshow(thmid_m2, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  #plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('mid-point')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(223)
  vmin0=-16
  vmax0=-12
  thneg_m1=np.ma.array(thneg,mask=peak<trmin)
  thneg_m2=np.ma.array(thneg_m1,mask=peak>trmax)
  imgplot = plt.imshow(thneg_m2, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-16,-15,-14,-13,-12],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('left slope')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(224)
  vmin0=12
  vmax0=16
  thpos_m1=np.ma.array(thpos,mask=peak<trmin)
  thpos_m2=np.ma.array(thpos_m1,mask=peak>trmax)
  imgplot = plt.imshow(thpos_m2, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[12,13,14,15,16],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('right slope')
  plt.xticks(xti)
  plt.yticks(yti)




     