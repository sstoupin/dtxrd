from pylab import *
import numpy as np
r2d = 180.0/np.pi
def figplot(s,dx,dy,indx1,indx2,indy1,indy2,peak,fwhm,stdev,com,thmid,thneg,thpos):
##
  thc=21.96
  thsi=23.6449
  dth = 2.0*(thsi-thc)
  dy = dy*cos(dth/r2d)/sin(thc/r2d)
  xyrange=(0.0,dx*(indx2-indx1),0.0,dy*(indy2-indy1))
                                  
  matplotlib.rcParams.update({'font.size': 20})
  matplotlib.rcParams['xtick.major.pad']='8'
  matplotlib.rcParams['ytick.major.pad']='8'
  
  f31=plt.figure(figsize=(16,12))
  shr1=0.8
  pad1=0.02
  asp1=20.0
  lbp1=20.0
  xti = [0,1,2,3,4,5,6,7,8,9,10,11]
  yti = [0,1,2,3,4,5,6,7,8,9,10,11]
  #--------------------------------------------------------------------------------------
  ysize=xyrange[3] #; print(ysize)
  xsize=xyrange[1] #; print(xsize)
  #print("y/x = ", ysize/xsize)
  #--------------------------------------------------------------------------------------
  # rejection threshold, rctopo assigns a high positive value 1e9 to 
  # topograph pixels that fall outside of the specified dynamic range 
  # note that rctopo-fast doesn't do this, the mask should be based on the 
  # user defined thresholds (e.g. based on low/high values of the  
  # 'peak' array which contains intensity topograph)
  trmax = 10.0
  #--------------------------------------------------------------------------------------
  plt.subplot(221)
  vmin0=26.5
  vmax0=29.5
  fwhm_m=np.ma.array(fwhm,mask=peak>trmax)
  imgplot = plt.imshow(fwhm_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
# imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[26.5,27.0,27.5,28.0,28.5,29.0,29.5],shrink=shr1,pad=pad1,aspect=asp1)
  #plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('FWHM')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(222)
  vmin0=-2.0
  vmax0=2.0
  thmid_m=np.ma.array(thmid,mask=peak>trmax)
  imgplot = plt.imshow(thmid_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  #plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('mid-point')
  plt.xticks(xti)
  plt.yticks(yti)
  #--------------------------------------------------------------------------------------- 
  plt.subplot(223)
  vmin0=-2.0
  vmax0=2.0
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
  plt.subplot(224)
  vmin0=-2.0
  vmax0=2.0
  thpos_m=np.ma.array(thpos,mask=peak>trmax)
  imgplot = plt.imshow(thpos_m, aspect='equal', extent=xyrange, vmin=vmin0, vmax=vmax0)
#  imgplot.set_cmap('jet')  
  plt.colorbar(imgplot,ticks=[-2,-1,0,1,2],shrink=shr1,pad=pad1,aspect=asp1)
  plt.xlabel('x [mm]',labelpad=lbp1)
  plt.ylabel('y [mm]',labelpad=lbp1)
  plt.title('right slope')
  plt.xticks(xti)
  plt.yticks(yti)
     