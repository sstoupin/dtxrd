#!/usr/bin/env python

'''
a subroutine to calculate parameters of a Bragg reflection

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

from numpy import *
from okada_si import *
from stoupin_c import *
from carr_ge import *
from lucht_sph import *
from fh import *

from dtxrd0 import *
from dtxrd_k import *
from constants import *
from curvestat import *
#----------------------------------------------------------------------
# Constants
#----------------------------------------------------------------------
#r2d=180.0/pi	     # radians to degrees conversion
#----------------------------------------------------------------------
#hpl=4.13566733e-15     # [eV*s] Planck constant
#hpl_bar=1.05457266e-34 # [J*s] Planck constant
#qe=1.60217653e-19     # [C] Electron charge
#me=9.0193897e-31     # [kg] mass of electron
#cl=299792458.0e10      # [Angstrom/s] speed of light 
#re=2.8179402894e-5  #[A] classical radius of electron
#----------------------------------------------------------------------


def Ebragg_c(element,h,k,l,T):
   
   h=float(h)
   k=float(k)
   l=float(l)
   T=float(T)
   
   if element=='Si':  a=a_okada(T)      
   elif element=='C': a=a_stoupin(T)
   elif element=='Ge': a=a_carr(T)
   else: fatalError('element is either C, Si or Ge at the moment')

   dh=a/sqrt(h**2.0+k**2.0+l**2.0)                               
   Eb=0.5*hpl*cl/dh
   
   qx=0.5/dh   
   f0h=f0h_ICD(element,qx)
   f00=f0h_ICD(element,0) # at qx=0
   expF0 = expFh_dia(0,0,0)
   expF  = expFh_dia(h,k,l)
   sigh=debye_sears(element,T,qx)
   V=a**3.0
      
   return [a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element]           

def Ebragg_sph(element,h,k,l,T):
    if element=='Al2O3':
          a=a_lucht(T); c=c_lucht(T)
    else: fatalError('element is Al2O3 only at the moment')
    
    dh=a*c/sqrt(4.0/3.0*c**2.0*(h**2.0+k**2.0+h*k)+a**2.0*l**2.0)
    Eb=0.5*hpl*cl/dh
    
    qx=0.5/dh
    f0h=1.0
    f00=1.0
    expF0 = 1.0
    expF = 1.0
    sigh = 1.0
    V=a**2.0*c
    
    return [[a,c],dh,Eb,f0h,f00,expF0,expF,sign,V,element]    
    

def thc_find(Ex,eta,phi,dc,crystal,P):               
   Eb=crystal[2]   
   thb=arcsin(Eb/Ex)   
   thc=thb; ee=1.0; i=0; imax=100
   #
   while (ee > 1.0e-6):
      i=i+1                       #; print i
      if i > imax: 
         print ('thc_find convergence problem!')
         break
      else:      
        k=[cos(thc)*cos(phi),cos(thc)*sin(phi),-sin(thc)]
        k_pr,thx,dth,eps,bh,Ty,Ry=dtxrd_k(k,eta,dc,Ex,P,crystal)        
        ee=abs(thx-thc)/dth      #; print eps
        thc=thx
        if 1.0-P < 1.0e-6:
           P=cos(2.0*thx)
        thc_pr=arcsin(k_pr[2]) 
        
   return [thc,thc_pr,eps,bh,Ty,Ry]

def thmax_find(Ex,eta,phi,dc,crystal,P,flag):               

   thc,thc_pr,eps,bh,Ty,Ry=thc_find(Ex,eta,phi,dc,crystal,P)      
#  Chi,wh_s,wh,eps_s,eps,dth_s,dth,de,Tplot,Rplot,dth_pr,eps_pr
   result0=dtxrd0(thc,eta,phi,dc,Ex,P,crystal)
   Chi=result0[0]
   dth=result0[6]
   thv=thc+1.5*dth*arange(-1.0,1.0,2.0/1000.0)
   
   T=[]; R=[]
   for thx in thv:   
     k=[cos(thx)*cos(phi),cos(thx)*sin(phi),-sin(thx)]
     [k_pr,Tx,Rx]=dtxrd1_k(k,eta,dc,Ex,P,crystal,Chi)
     T=T+[Tx]; R=R+[Rx]

   thv=array(thv); R=array(R); T=array(T)  
   if flag=='R':   stat1=curvestat(thv,R,0.0)
   elif flag=='T': stat1=curvestat(thv,T,0.0)
   else: stat1=0; print " flag input error "
   
   thmax=stat1[0]
   k=[cos(thmax)*cos(phi),cos(thmax)*sin(phi),-sin(thmax)]
   [k_pr,Tx,Rx]=dtxrd1_k(k,eta,dc,Ex,P,crystal,Chi)                       
   thmax_pr=arcsin(k_pr[2])      
   
   return [thmax,thmax_pr,eps,bh,Ty,Ry]
               
def thr_find(Ex,eta,phi,dc,crystal,P):
   a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element=crystal
   E0=Ex; ee=1.0; i=0; imax=100
   #
   while (ee > 1.0e-16):
      i=i+1
      if i > imax:
         print ('thr_find convergence problem!')
         break
      else:   
         fa=fa_asf(element,E0)
         F0 = (f00+fa)*expF0
         #Chi0  = -re*F0/(pi*V)*(hpl*cl/E0)**2.0
         #wh_s  = -2.0*real(Chi0)*(dh/lamx)**2.0    
         wh_s=-2.0*real(-re*F0/(pi*V))*dh**2.0
         #print 'wh_s = ', wh_s
         # theta_r:
         denr=1.0+sqrt(1.0+4.0*wh_s*(tan(eta))**2.0)   
         bthr=arctan(2.0*wh_s*tan(eta)/denr)   
         thr=0.5*pi-bthr
         # lambda_r:
         lamr=2.0*dh*cos(bthr)*(1.0-tan(bthr)/tan(eta))
         Er=hpl*cl/lamr
         ee=abs(E0-Er)/Er      
         E0=Er
   k=[cos(thr)*cos(phi),cos(thr)*sin(phi),-sin(thr)]
   k_pr,thx,dth,eps,bh,Ty,Ry=dtxrd_k(k,eta,dc,Er,P,crystal)        
   thr_pr=arcsin(k_pr[2]) 
               
   return [thr,Er,eps,bh,Ty,Ry]
   
