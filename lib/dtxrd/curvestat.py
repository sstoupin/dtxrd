'''
subroutines to calculate statistical parameters of a curve (FWHM, COM, STDEV, etc.)

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

from numpy import *
#from scipy.optimize import *
#from pylab import *
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline

###################################################################################
# Basic Shapes
###################################################################################

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
    return b0 + (b[1]-b0)/(1+(x-b[2])**2.0/b[3]**2.0)    
    #return b[0] + (b[1]-b[0])*b[3]/((x-b[2])**2.0 + b[3]**2.0) # this is canonical lorentzian but it doesn't matter
    # once we define the amplitude without renormalization

# without background: 
def gauss0(a,x):
    # here the input a[3] = sqrt(2.0)*sigma 
    # usually fwhm = 2.0*sqrt(2.0 * log(2.0))*sigma
    # a[3] = sqrt(2.0) * fwhm/(2.0*sqrt(2.0 * log(2.0)) = 0.5*fwhm/(sqrt(log(2.0))
    ###########################
    #a0=abs(a[0])
    return a[0]*exp(-(x-a[1])**2.0/a[2]**2.0)

# without background:
def lorentz0(b,x):
    # fwhm = 2*b[2] for canonical lorentzian
    #b0=abs(b[0])
    return b[0]/(1+(x-b[1])**2.0/b[2]**2.0)
    #return b[0]*b[2]/((x-b[1])**2.0 + b[2]**2.0)    # this is canonical lorentzian
#######################################################################################
# Point-by-point methods
#######################################################################################

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


def curvestat0(th,r,bkg):  # expedited procedure, no interpolation - only approx. answer    
    r_max=max(r)
    th_max = th[(r == r_max)]
    th_max = th_max[0]
    half=0.5*(r_max-bkg)+bkg
    # 
    inds = [x for x in range(len(r)) if r[x] > half]
    th_neg = th[min(inds)]
    th_pos = th[max(inds)]
    #                        
    th_mid=0.5*(th_neg+th_pos)
    fwhm=abs(th_pos-th_neg)
    #
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
    int=norm*(th[len(th)-1] - th[0])/float((len(th)-1))
    #int = norm
    #
    return [th_max,r_max,th_neg,th_pos,th_mid,fwhm,com,var,int]


#################################################################################################
## Vectorized methods
#################################################################################################

def curvestat_(th,rcurve):
    '''
    Find curve statistics using zero crossings 
    where data is 3D array sorted such that
    data[0] is the angular dependent value
        
    on 03/28/2010 replaced "gradient" method with diff (consecutive difference) for interpolation
    also in finding thneg revised to start from the -1 point below the zero-crossing
    this yields results similar to rctopo point-by-point algorithm
    '''
    rmax_ = amax(rcurve, axis=0)
    const = ones(rcurve.shape[0],dtype=float32)
    th_ = th[...,newaxis,newaxis]*ones(rcurve.shape)
    rcurve_eff = where(rcurve<0.0, 0.0, rcurve)  # reject negative values!
    #th_eff = th_
    norm = sum(abs(rcurve_eff),axis=0)  # 2D array
    com = sum(th_*rcurve_eff,axis=0)/norm # need to exclude rcurve_ - bkg <0        
    com_ = ones(rcurve_eff.shape)*com[newaxis,...]  # extend to 3D
    var_ = abs(rcurve_eff)*(th_ - com_)**2.0
    stdev = sqrt(sum(var_, axis=0)/norm)
    # 09/17/2019
    # success!
    # to do: regrid index to ugol scale in vector form
    # interpolation possible?
    # 09/23/2019 interpolation success
    # why inflated rms and pv?
    #------------------------------------------------------------------------------------------------
    # subtract 0.5 of max to form a curve (3D) where zero crossings correspond to left/right slopes
    dfwhm_ = rcurve - 0.5*const[...,newaxis,newaxis]*rmax_[newaxis,...]    
    #print('dfwhm_ = ', dfwhm_)
    dfwhm_m = ma.masked_where(dfwhm_ < 0.0,dfwhm_) # mask negative values
    #fwhm = ma.count(dfwhm_m,axis=0) # could be but not
    #-----------------------
    # derivative
    #dr_ = gradient(dfwhm_,axis=0)  # 3D derivative/differences for intensity
    dr_ = diff(dfwhm_,axis=0,append=0.0) #dfwhm_[-1,:,:])
    #dth_ = gradient(th_,axis=0)    # 3D derivative/differences for the angular scale  
    dth_ = diff(th_,axis=0,append=1.0)   #th_[-1,:,:]) # th is in microradians at this point
    derr_ = dr_/dth_               # 3D the actual derivative
    # arrays of indicies at left/right edges as notmasked_edges
    thneg0,thpos0 = asarray(ma.notmasked_edges(dfwhm_m,axis=0))
    # initialize arrays:
    thneg2 = zeros((rcurve.shape[1],rcurve.shape[2])) 
    thpos1 = zeros((rcurve.shape[1],rcurve.shape[2]))        
    iz2 = zeros((rcurve.shape[1],rcurve.shape[2]),dtype=int32)
    iz1 = zeros((rcurve.shape[1],rcurve.shape[2]),dtype=int32)
    #
    z1,x1,y1 = thneg0[0],thneg0[1],thneg0[2]    
    thneg2[x1,y1] = take(th,z1-1)             # z1 here are indicies 
    iz1[x1,y1]  = z1-1                          #2D array of indicies
    z2,x2,y2 = thpos0[0],thpos0[1],thpos0[2]    
    thpos1[x2,y2] = take(th,z2)             #z2 here are indicies
    iz2[x2,y2] = z2                           #2D array of indicies
    #
    #print('iz1 = ', iz1)
    #print('iz2 = ', iz2)
    # -------------------------------------------------------------------
    # extracting 2D arrays from 3D based on 2D index arrays iz1 and iz2
    m,n = iz1.shape
    I,J = ogrid[:m,:n]    
    rneg2 = dfwhm_[iz1,I,J]    
    rpos1 = dfwhm_[iz2,I,J]
    derneg2 = derr_[iz1,I,J]
    derpos1 = derr_[iz2,I,J]
    #---------------------------------------------------------------------
    #
    eps0 = 1.0e-3 # ct/urad    # this is to avoid division by zero
    thneg = thneg2 - rneg2/(derneg2 + eps0)
    thpos = thpos1 - rpos1/(derpos1 + eps0)
    ### these are derived metrics          
    fwhm = thpos-thneg
    thmid = 0.5*(thpos + thneg)  
    #
    return [thneg,thpos,thmid,fwhm,com,stdev]