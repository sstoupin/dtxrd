#!/usr/bin/env python

'''
a subroutine to calculate susceptibility(Chi) for a Bragg reflection

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''
###########################################################################################
# v 0.03 #######################################################
###########################################################################################
from numpy import *
#from okada_si import *
#from stoupin_c import *
#from carr_ge import *
#from lucht_sph import *
from fh import *
#
#from dtxrd import *
#from dtxrd_k import *
from constants import *
#from curvestat import *
import sys
######################################################################################
def fatalError(msg):
  sys.stderr.write('Error: ')
  sys.stderr.write(str(msg))
  sys.stderr.write('\n')
  sys.exit(1)                                

def chi(element,h,k,l,T,Ex):
  epsFh = 1.0e-6
  flagFh = 1

  h = float(h)
  k = float(k)
  l = float(l)
  T = float(T)
  Ex = float(Ex)
  lamx = hpl*cl/Ex
  
  if element=='Si':  
     from okada_si import a_okada
     a=a_okada(T)
     #    
     dh=a/sqrt(h**2.0+k**2.0+l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh
     #
     f0h=f0h_ICD(element,qx)
     f00=f0h_ICD(element,0.0) # at qx=0
     expF0 = expFh_dia(0,0,0)
     expF  = expFh_dia(h,k,l) #; print "expF(h,k,l) = ", expF
     B_el, sigh = debye_sears(element,T,qx)
     V=a**3.0
     #
     fa=fa_asf(element,Ex)
     Fh = sigh*(f0h+fa)*expF                  #; print Fh
     Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
     F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???
     #
     if abs(expF) < epsFh: flagFh = 0
  #----------------------------------------------------------------------------   
  elif element=='C': 
     from stoupin_c import a_stoupin
     a=a_stoupin(T)
     #    
     dh=a/sqrt(h**2.0+k**2.0+l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh
     #
     f0h=f0h_ICD(element,qx)
     f00=f0h_ICD(element,0.0) # at qx=0
     expF0 = expFh_dia(0,0,0)
     expF  = expFh_dia(h,k,l) #; print "expF(h,k,l) = ", expF
     B_el, sigh = debye_sears(element,T,qx)
     V=a**3.0
     #
     fa=fa_asf(element,Ex)
     Fh = sigh*(f0h+fa)*expF                  #; print Fh
     Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
     F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???
     #
     if abs(expF) < epsFh: flagFh = 0
  #-----------------------------------------------------------------------------    
  elif element=='Ge': 
     from carr_ge import a_carr
     a=a_carr(T)
     #    
     dh=a/sqrt(h**2.0+k**2.0+l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh
     #
     f0h=f0h_ICD(element,qx)
     f00=f0h_ICD(element,0.0) # at qx=0
     expF0 = expFh_dia(0,0,0)
     expF  = expFh_dia(h,k,l) #; print "expF(h,k,l) = ", expF
     B_el, sigh = debye_sears(element,T,qx)
     V=a**3.0
     #
     fa=fa_asf(element,Ex)
     Fh = sigh*(f0h+fa)*expF                  #; print Fh
     Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
     F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???
     #
     if abs(expF) < epsFh: flagFh = 0
  #------------------------------------------------------------------------------                             
  elif element=='Al2O3':           
     from lucht_sph import a_lucht, c_lucht
     a=a_lucht(T); c=c_lucht(T)               
  
     dh=a*c/sqrt(4.0/3.0*c**2.0*(h**2.0+k**2.0+h*k)+a**2.0*l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh    
     
     f0h_Al=f0h_ICD('Al',qx); f0h_O=f0h_ICD('O',qx)
     f00_Al=f0h_ICD('Al',0.0); f00_O=f0h_ICD('O',0.0)
     expF0_Al = expFh_sph('Al',0,0,0); expF0_O = expFh_sph('O',0,0,0); 
     expF_Al = expFh_sph('Al',h,k,l) #; print "expF_Al(h,k,l) =  ", expF_Al
     expF_O = expFh_sph('O',h,k,l)   #; print "expF_O(h,k,l)  =  ", expF_O
     
     if abs(expF_Al)<epsFh and abs(expF_O)<epsFh: flagFh=0
     
     B_Al=0.195 # Angstrom^2    !!! temp dependece?
     B_O=0.274  # Angstrom^2    !!! temp dependence?
     sigh_Al=exp(-B_Al*qx**2.0)
     sigh_O=exp(-B_O*qx**2.0)
     V=sqrt(3.0)/2.0*a**2.0*c
     #
     fa_Al=fa_asf('Al',Ex); fa_O=fa_asf('O',Ex)
     Fh = sigh_Al*(f0h_Al+fa_Al)*expF_Al+sigh_O*(f0h_O+fa_O)*expF_O
     Fh_= sigh_Al*(f0h_Al+fa_Al)*expF_Al.conjugate()+sigh_O*(f0h_O+fa_O)*expF_O.conjugate()
     F0 = (f00_Al+fa_Al)*expF0_Al+(f00_O+fa_O)*expF0_O
  #------------------------------------------------------------------------------------------
  elif element=='FeS2':
     from chrystall_pyr import a_chrystall
     a=a_chrystall(T)
     #
     dh=a/sqrt(h**2.0+k**2.0+l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh    
     #
     f0h_Fe=f0h_ICD('Fe',qx); f0h_S=f0h_ICD('S',qx)
     f00_Fe=f0h_ICD('Fe',0.0); f00_S=f0h_ICD('S',0.0)
     expF0_Fe = expFh_pyr('Fe',0,0,0); expF0_S = expFh_pyr('S',0,0,0); 
     expF_Fe = expFh_pyr('Fe',h,k,l) #; print "expF_Al(h,k,l) =  ", expF_Al
     expF_S = expFh_pyr('S',h,k,l)   #; print "expF_O(h,k,l)  =  ", expF_O
     #
     if abs(expF_Fe)<epsFh and abs(expF_S)<epsFh: flagFh=0
     #
     B_Fe = 0.28 # Angstrom^2    !!! temp dependece?   from Bayliss_AM1977
     B_S = 0.33  # Angstrom^2    !!! temp dependence?  from Bayliss_AM1977
     sigh_Fe = exp(-B_Fe*qx**2.0)
     sigh_S = exp(-B_S*qx**2.0)          
     V = a**3.0          
     #
     fa_Fe=fa_asf('Fe',Ex); fa_S=fa_asf('S',Ex)     
     Fh = sigh_Fe*(f0h_Fe+fa_Fe)*expF_Fe+sigh_S*(f0h_S+fa_S)*expF_S     
     Fh_= sigh_Fe*(f0h_Fe+fa_Fe)*expF_Fe.conjugate()+sigh_S*(f0h_S+fa_S)*expF_S.conjugate()     
     F0 = (f00_Fe+fa_Fe)*expF0_Fe+(f00_S+fa_S)*expF0_S
  #-------------------------------------------------------------------------------------------   
  elif element=='SiC-4H':           
     from springer_sic import a_SiC4H, c_SiC4H
     a=a_SiC4H(T); c=c_SiC4H(T)               
  
     dh=a*c/sqrt(4.0/3.0*c**2.0*(h**2.0+k**2.0+h*k)+a**2.0*l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh    
     
     f0h_Si=f0h_ICD('Si',qx); f0h_C=f0h_ICD('C',qx)
     f00_Si=f0h_ICD('Si',0.0); f00_C=f0h_ICD('C',0.0)
     expF0_Si = expFh_SiC4H('Si',0,0,0); expF0_C = expFh_SiC4H('C',0,0,0); 
     expF_Si = expFh_SiC4H('Si',h,k,l) #; print "expF_Al(h,k,l) =  ", expF_Al
     expF_C = expFh_SiC4H('C',h,k,l)   #; print "expF_O(h,k,l)  =  ", expF_O
     
     if abs(expF_Si)<epsFh and abs(expF_C)<epsFh: flagFh=0
     
     # these are from SiC-6H (Capitani) - need to find better approximation     
     ### Anisotropic displacement parameters from Capitani etal AM 2007
     U11=5.0e-3 # Angstrom^-2
     U22=5.0e-3 # Angstrom^-2
     U33=5.0e-3 # Angstrom^-2     
     U12=2.0e-3 # Angstrom^-2
     e11=U11*h**2.0*4.0/(3.0*a**2.0)  # reciprocal lattice constants sqrt(4)/(sqrt(3)*a)
     e22=U22*k**2.0*4.0/(3.0*a**2.0)  # sqrt(4)/(sqrt(3)*a)
     e33=U33*l**2.0/c**2.0            # 1/c
     e12=U12*h*k*4.0/(3.0*a**2.0)
     sigh_Si=exp(-2.0*pi**2.0*(e11+e22+e33+2.0*e12))
     sigh_C=sigh_Si               
     V=sqrt(3.0)/2.0*a**2.0*c
     #
     fa_Si=fa_asf('Si',Ex); fa_C=fa_asf('C',Ex)
     Fh = sigh_Si*(f0h_Si+fa_Si)*expF_Si+sigh_C*(f0h_C+fa_C)*expF_C
     Fh_= sigh_Si*(f0h_Si+fa_Si)*expF_Si.conjugate()+sigh_C*(f0h_C+fa_C)*expF_C.conjugate()
     F0 = (f00_Si+fa_Si)*expF0_Si+(f00_C+fa_C)*expF0_C
  #-------------------------------------------------------------------------------------------   
  elif element=='SiC-6H':           
     from springer_sic import a_SiC6H, c_SiC6H
     a=a_SiC6H(T); c=c_SiC6H(T)               
  
     dh=a*c/sqrt(4.0/3.0*c**2.0*(h**2.0+k**2.0+h*k)+a**2.0*l**2.0)
     Eb=0.5*hpl*cl/dh
     qx=0.5/dh    
     
     f0h_Si=f0h_ICD('Si',qx); f0h_C=f0h_ICD('C',qx)
     f00_Si=f0h_ICD('Si',0.0); f00_C=f0h_ICD('C',0.0)
     expF0_Si = expFh_SiC6H('Si',0,0,0); expF0_C = expFh_SiC6H('C',0,0,0); 
     expF_Si = expFh_SiC6H('Si',h,k,l) #; print "expF_Al(h,k,l) =  ", expF_Al
     expF_C = expFh_SiC6H('C',h,k,l)   #; print "expF_O(h,k,l)  =  ", expF_O
     
     if abs(expF_Si)<epsFh and abs(expF_C)<epsFh: flagFh=0
     
     #B_Si=0.195 # Angstrom^2    !!! temp dependece?
     #B_C=0.274  # Angstrom^2    !!! temp dependence?     
     #sigh_Si=exp(-B_Si*qx**2.0)
     #sigh_C=exp(-B_C*qx**2.0)
     ### Anisotropic displacement parameters from Capitani etal AM 2007
     U11=5.0e-3 # Angstrom^-2
     U22=5.0e-3 # Angstrom^-2
     U33=5.0e-3 # Angstrom^-2     
     U12=2.0e-3 # Angstrom^-2
     e11=U11*h**2.0*4.0/(3.0*a**2.0)  # reciprocal lattice constants sqrt(4)/(sqrt(3)*a)
     e22=U22*k**2.0*4.0/(3.0*a**2.0)  # sqrt(4)/(sqrt(3)*a)
     e33=U33*l**2.0/c**2.0            # 1/c
     e12=U12*h*k*4.0/(3.0*a**2.0)
     sigh_Si=exp(-2.0*pi**2.0*(e11+e22+e33+2.0*e12))
     sigh_C=sigh_Si     
     V=sqrt(3.0)/2.0*a**2.0*c
     #
     fa_Si=fa_asf('Si',Ex); fa_C=fa_asf('C',Ex)
     Fh = sigh_Si*(f0h_Si+fa_Si)*expF_Si+sigh_C*(f0h_C+fa_C)*expF_C
     Fh_= sigh_Si*(f0h_Si+fa_Si)*expF_Si.conjugate()+sigh_C*(f0h_C+fa_C)*expF_C.conjugate()
     F0 = (f00_Si+fa_Si)*expF0_Si+(f00_C+fa_C)*expF0_C
  #-------------------------------------------------------------------------------------------     
  else:
     fatalError('available elements are C, Si, Ge, Al2O3, FeS2, SiC-4H and SiC-6H')
  #     
  Chi0  = -re*F0/(pi*V)*lamx**2.0   #; print "chi_{0} = ",  Chi0
  Chih  = -re*Fh/(pi*V)*lamx**2.0   #; print "chi_{h} = ",  Chih
  Chih_ = -re*Fh_/(pi*V)*lamx**2.0  #; print "chi_{-h} = ", Chih_
                                         
  Chi=[Chi0,Chih,Chih_]   
  
  return [[Chi,dh],flagFh]
    