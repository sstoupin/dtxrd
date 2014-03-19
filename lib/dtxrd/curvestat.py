'''
a subroutine to calculate statistical parameters for a reflectivity curve (FWHM, COM, etc.)

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

from numpy import *
#from scipy.optimize import *
#from pylab import *

def curvestat(th,r,bkg):
      # th - argument
      # r - data
      # bkg -background
      
      N1=len(th); N2=len(r)
      if N1 !=N2:
        print 'Error: length of data vector does not match that of the argument!'
      else: 
        th_l=list(th); r_l=list(r)
        i_max=r_l.index(max(r))                                              
        th_max=th[i_max]
        r_max=max(r)

        half=0.5*(r_max-bkg)+bkg
        x0=r[0]
        for x in r:                
                if x>half:
                    if x0<half:
                      i_1=r_l.index(x0); i_2=r_l.index(x)
                      r1=r[i_1]; r2=r[i_2]
                      th1=th[i_1]; th2=th[i_2]
                      der1=(r2-r1)/(th2-th1) # linear interpolation: r=der1*(th-th1)+r1                      
                      th_neg=th1+(half-r1)/der1
                      
                elif x<half:
                    if x0>half:
                      i_1=r_l.index(x0); i_2=r_l.index(x)
                      r1=r[i_1]; r2=r[i_2]
                      th1=th[i_1]; th2=th[i_2]
                      der2=(r2-r1)/(th2-th1) # linear interpolation: r=der2*(th-th1)+r1        
                      th_pos=th1+(half-r1)/der2
                x0=x                                                       
                        
        th_mid=0.5*(th_neg+th_pos)
        fwhm=abs(th_pos-th_neg)
        com=sum(th*(r-bkg))/sum(r-bkg)
        return [th_max,r_max,th_neg,th_pos,th_mid,fwhm,com]

def gauss(a,x):
    a0=abs(a[0])
    return a0+(a[1]-a0)*exp(-(x-a[2])**2.0/a[3]**2.0)

def lorentz(b,x):
    b0=abs(b[0])
    return b0+(b[1]-b0)/(1+(x-b[2])**2.0/b[3]**2.0)
