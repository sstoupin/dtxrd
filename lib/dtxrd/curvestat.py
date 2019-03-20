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
from scipy.integrate import simps

def curvestat(th,r,bkg):
      # th - argument
      # r - data
      # bkg -background
      
      N1=len(th); N2=len(r)
      if N1 !=N2:
        print('Error: length of data vector does not match that of the argument!')
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
        
        r_eff0 = r - bkg
        r_eff = r_eff0[r_eff0 > 0] 
        th_eff = th[r_eff0 > 0]
        #        
        #norm = simps(r_eff,th_eff)
        #!!!!!!!!!!!! SIMPS failed on 23OCT2017 - sample NDT111-3_Bc2!!!!!!!        
        norm = sum(r_eff)
        com=sum(th_eff*(r_eff))/norm             # center of mass (mean) or  first cumulant
        #com = simps(th_eff*r_eff,th_eff)/norm
        #var=sum((r_eff)*(th_eff-com)**2.0)/sum(r_eff)  # variance or second cumulant
        var_num = r_eff*(th_eff-com)**2.0        
        #var = simps(var_num,th_eff)/norm        
        #var = simps(th_eff**2.0*r_eff,th_eff)/norm - com**2.0                
        #var=abs(sum((r-bkg)*th**2.0)/sum(r-bkg)-com**2.0)  # equivalent!
        var=sum(var_num)/norm                
        int=norm*(th[len(th)-1] - th[0])/float((N1-1))
        #int = norm
        #
        return [th_max,r_max,th_neg,th_pos,th_mid,fwhm,com,var,int]

def gauss(a,x):
    # here the input a[3] = sqrt(2.0)*sigma 
    # usually fwhm = 2.0*sqrt(2.0 * log(2.0))*sigma
    # a[3] = sqrt(2.0) * fwhm/(2.0*sqrt(2.0 * log(2.0)) = 0.5*fwhm/(sqrt(log(2.0))
    ###########################
    a0=abs(a[0])
    return a0+(a[1]-a0)*exp(-(x-a[2])**2.0/a[3]**2.0)

def lorentz(b,x):
    # fwhm = 2*b[3]
    b0=abs(b[0])
    return b0+(b[1]-b0)/(1+(x-b[2])**2.0/b[3]**2.0)
