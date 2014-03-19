#!/usr/bin/python

################################################################
## subroutine to calculate susceptibility of hkl reflection
# v 0.02 #######################################################
# Stanislav Stoupin ## sstoupin@aps.anl.gov ####################
################################################################

#from numpy import *
from okada_si import *
from stoupin_c import *
from carr_ge import *
from lucht_sph import *
from fh import *
#from dtxrd import *
#from dtxrd_k import *
from constants import *
#from curvestat import *
######################################################################################
def chi(element,h,k,l,T,Ex):
  epsFh=1.0e-6
  flagFh=1

  h=float(h)
  k=float(k)
  l=float(l)
  T=float(T)
  Ex=float(Ex)
  
  if element=='Si':  
    a=a_okada(T)
  elif element=='C': 
    a=a_stoupin(T)
  elif element=='Ge': 
    a=a_carr(T)
  elif element=='Al2O3':
    a=a_lucht(T); c=c_lucht(T)               
  else: fatalError('element is either C, Si, Ge or Al2O3 at the moment')

  lamx=hpl*cl/Ex
      
  if element=='Si' or element=='C' or element=='Ge':
   
     dh=a/sqrt(h**2.0+k**2.0+l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh
     
     f0h=f0h_ICD(element,qx)
     f00=f0h_ICD(element,0.0) # at qx=0
     expF0 = expFh_dia(0,0,0)
     expF  = expFh_dia(h,k,l) #; print "expF(h,k,l) = ", expF
     sigh=debye_sears(element,T,qx)
     V=a**3.0
     
     fa=fa_asf(element,Ex)
     Fh = sigh*(f0h+fa)*expF                  #; print Fh
     Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
     F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???
     
     if abs(expF) < epsFh: flagFh = 0
                             
  elif element=='Al2O3':
     dh=a*c/sqrt(4.0/3.0*c**2.0*(h**2.0+k**2.0+h*k)+a**2.0*l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh    
     
     f0h_Al=f0h_ICD('Al',qx); f0h_O=f0h_ICD('O',qx)
     f00_Al=f0h_ICD('Al',0.0); f00_O=f0h_ICD('O',0.0)
     expF0_Al = expFh_sph('Al',0,0,0); expF0_O = expFh_sph('O',0,0,0); 
     expF_Al = expFh_sph('Al',h,k,l) #; print "expF_Al(h,k,l) =  ", expF_Al
     expF_O = expFh_sph('O',h,k,l)   #; print "expF_O(h,k,l)  =  ", expF_O
     
     if abs(expF_Al)<epsFh and abs(expF_O)<epsFh: flagFh=0
     
     B_Al=0.195 # Angstrom^2
     B_O=0.274  # Angstrom^2
     sigh_Al=exp(-B_Al*qx**2.0)
     sigh_O=exp(-B_O*qx**2.0)
     V=sqrt(3.0)/2.0*a**2.0*c
     fa_Al=fa_asf('Al',Ex); fa_O=fa_asf('O',Ex)
     Fh = sigh_Al*(f0h_Al+fa_Al)*expF_Al+sigh_O*(f0h_O+fa_O)*expF_O
     Fh_= sigh_Al*(f0h_Al+fa_Al)*expF_Al.conjugate()+sigh_O*(f0h_O+fa_O)*expF_O.conjugate()
     F0 = (f00_Al+fa_Al)*expF0_Al+(f00_O+fa_O)*expF0_O
     
     
  Chi0  = -re*F0/(pi*V)*lamx**2.0   #; print "chi_{0} = ",  Chi0
  Chih  = -re*Fh/(pi*V)*lamx**2.0   #; print "chi_{h} = ",  Chih
  Chih_ = -re*Fh_/(pi*V)*lamx**2.0  #; print "chi_{-h} = ", Chih_
                                         
  Chi=[Chi0,Chih,Chih_]   
  
  return [[Chi,dh],flagFh]
    
