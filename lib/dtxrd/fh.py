#!/usr/bin/env python
########################################################
## SET OF FUNCTIONS TO CALCULATE STRUCTURE FACTORS
########################################################
# VER=0.18
#########################################################################
# to do list:
# 1. check alternatives for DWF
#########################################################################

#GLOBAL:
import dtxrd
import os
libpath = os.path.dirname(dtxrd.__file__)

#LOCAL:
#import sys
#import commands
#cmdout=commands.getstatusoutput('echo $HOME')
#libpath=cmdout[1]+'/bin/DTXRD'
#sys.path.append(libpath)

from numpy import *
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import KroghInterpolator
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import interp1d

hpl=4.13566733e-15 # [eV*s] Planck constant
c=299792458e+10    # [Angstrom/s] speed of light

###############################################################
### ATOMIC FORM FACTOR
###############################################################

def f0h_ICD(atom,qx):
#    lamx=hpl*c/Ex   #; print "lamx [A] = ", lamx
#    qx=sin(th)/lamx #; print "qx [A-1] = ", qx
        
    fileName=libpath+'/data/f0h_ICD.dat'
    data=open(fileName, 'r')
    while 1:
        line=data.readline()
        stuff=line.split(' ')
        if stuff[0]=='Q':
          line1=data.readline()
          q=line1.split(' ')
        elif stuff[0]==atom:
          line1=data.readline()
          f0=line1.split(' ')
          break
        elif not(line): 
          break
    
    q=array(q, dtype=float)   #; print q
    f0=array(f0, dtype=float) #; print f0
    f0q=interp1d(q,f0)  # UnivariateSpline KroghInterpolator BarycentricInterpolator  do not work here !!!!
                        # interp1d works fine but limited range!
    # TEST
    #figure(1)
    #plot(q,f0,'bo')
##    qplot=arange(0,6,0.01)
##    plot(qplot,f0q(qplot),'r-')
##    show()
    return f0q(qx)
#------------------------------------------------------------------------

def f0h_waasmaier(atom,qx):
#    lamx=hpl*c/Ex   #; print "lamx [A] = ", lamx
#    qx=sin(th)/lamx #; print "qx [A-1] = ", qx

    a=zeros(6); b=zeros(6)
    fileName=libpath+'/data/waasmaier.dat'
    data=open(fileName, 'r')
    while 1:
        line=data.readline()
        stuff=line.split(' ')                
        if stuff[0]==atom:
          a[1]=stuff[2];  b[1]=stuff[3]
          a[2]=stuff[4];  b[2]=stuff[5]
          a[3]=stuff[6];  b[3]=stuff[7]
          a[4]=stuff[8];  b[4]=stuff[9]
          a[5]=stuff[10]; b[5]=stuff[11]
          c=float(stuff[12])
          break
        elif not(line): 
          break
    a=array(a, dtype=float)
    b=array(b, dtype=float)
    
    f0q=[]
    for q in qx:
      f0q=f0q+[sum(a*exp(-b*q**2.0))+c]  
    print len(qx)
    print len(f0q)
    return array(f0q)

## TEST
## test results (03/17/2012):
## interp1d gives better fit to the points in the ICD tables
## the a5 b5 c sum of Waasmaier is a useful analytical form for the approximation
## it does depart a bit from the ICD tabulated data but possibly gives better than interp1d approximation in between the points
## for Si the difference between the two does not exceed 0.021
## for C  the difference between the two does not exceed 0.012
#from pylab import *
#qx=arange(0.0,6.0,0.01)
#f0h_icd=f0h_ICD('C',qx)
#f0h_waas=f0h_waasmaier('C',qx)
#figure(1)
#plot(qx,f0h_icd,'b-')
#plot(qx,f0h_waas,'r-')
#figure(2)
#diff=f0h_icd-f0h_waas
#plot(qx,diff)
#show()
#
############################################################################################################
###  ANOMALOUS CORRECTIONS
############################################################################################################
# dictionary format 
# f=elements{'atom'}
# Z=f[0]
# f_p=f[1]  Cromer_JPC1970 at wavelengths
# f_pp=f[2] Cromer_JPC1970 at wavelengths
# df=f[3]   Kissel_AC1990

wl=array([2.28962,1.93597,1.54052,0.70926,0.55936])
E=hpl*c/wl

elements = {
      'C'  : [6, [0.035,0.026,0.017,0.002,0.000],[0.021,0.015,0.009,0.002,0.001],0.001],
      'Si' : [14,[0.355,0.311,0.244,0.072,0.042],[0.693,0.509,0.330,0.071,0.043],0.011]
      }

elementsZ = {
	'H' : 1,
	'He': 2,
	'Li': 3,
	'Be': 4,
	'B' : 5,
	'C' : 6,
	'N' : 7,
	'O' : 8,
	'F' : 9,
	'Ne': 10,
	'Na': 11,
	'Mg': 12,
	'Al': 13,
	'Si': 14,
	'P' : 15,
	'S' : 16,
	'Cl': 17,
	'Ar': 18,
	'K' : 19,
	'Ca': 20,
	'Sc': 21,
	'Ti': 22,
	'V' : 23,
	'Cr': 24,
	'Mn': 25,
	'Fe': 26,
	'Co': 27,
	'Ni': 28,
	'Cu': 29,
	'Zn': 30,
	'Ga': 31,
	'Ge': 32,
	'As': 33,
	'Se': 34,
	'Br': 35,
	'Kr': 36,
	'Rb': 37,
	'Sr': 38,
	'Y' : 39,
	'Zr': 40,
	'Nb': 41,
	'Mo': 42,
	'Tc': 43,
	'Ru': 44,
	'Rh': 45,
	'Pd': 46,
	'Ag': 47,
	'Cd': 48,
	'In': 49,
	'Sn': 50,
	'Sb': 51,
	'Te': 52,
	'I' : 53,
	'Xe': 54,
	'Cs': 55,
	'Ba': 56,
	'La': 57,
	'Ce': 58,
	'Pr': 59,
	'Nd': 60,
	'Pm': 61,
	'Sm': 62,
	'Eu': 63,
	'Gd': 64,
	'Tb': 65,
	'Dy': 66,
	'Ho': 67,
	'Er': 68,
	'Tm': 69,
	'Yb': 70,
	'Lu': 71,
	'Hf': 72,
	'Ta': 73,
	'W' : 74,
	'Re': 75,
	'Os': 76,
	'Ir': 77,
	'Pt': 78,
	'Au': 79,
	'Hg': 80,
	'Tl': 81,
	'Pb': 82,
	'Bi': 83,
	'Po': 84,
	'At': 85,
	'Rn': 86,
	'Fr': 87,
	'Ra': 88,
	'Ac': 89,
	'Th': 90,
	'Pa': 91,
	'U' : 92,
	'Np': 93,
	'Pu': 94,
	'Am': 95,
	'Cm': 96,
	'Bk': 97,
	'Cf': 98,
	'Es': 99,
	'Fm': 100,
	'Md': 101,
	'No': 102,
	'Lr': 103,
	'Rf': 104,
	'Ha': 105,
	'Sg': 106,
	'Bh': 107,
	'Hs': 108,
	'Mt': 109,
	    }

def fa_cromer(atom,Ex):

    f=elements[atom]    
#    Efake=arange(E[0],E[4],0.01)
#    f_ppf=interp1d(E,array(f[2]))
#    f_ppfake=f_ppf(Efake)
#    f_p=UnivariateSpline(E,array(f[1]))        
#    f_pp=UnivariateSpline(E,array(f[2]))
    f_p=interp1d(E,array(f[1]))        
    f_pp=interp1d(E,array(f[2]))
    f_a=float(f[3])        
#    print "f_p = ", f_p(Ex)
#    print "f_pp = ", f_pp(Ex)  # things make sense only if taken with "-" sign - but not clear from the literature (Cromer, etc.)
          
    return f_p(Ex)+f_a-1.0j*f_pp(Ex)        
# --------------------------------------------------------------------------------    
def fa_asf(atom,Ex):
    from myio import readFile
    fileName=libpath+'/asf/'+atom+'.asf'
    d1,d2=readFile(fileName)
    Easf=d2[:,0]; Easf=1.0e3*Easf
    f1=d2[:,1]
    f2=d2[:,2]
    #
    Z=float(elementsZ[atom])
    f_1 = interp1d(Easf,f1)  # it seems that f1 in fa_asf is f1=f_p+f_rel+Z    
    f_2= interp1d(Easf,f2)
    #    
    return f_1(Ex)-Z+1.0j*f_2(Ex)
    
# Tests: 
##from pylab import * 
#from scipy.interpolate import UnivariateSpline 
#from scipy.interpolate import KroghInterpolator
#from scipy.interpolate import BarycentricInterpolator
##f=elements['Si']
##figure(1)
#plot(E,f[2],'ro')
##Ex=arange(5000,30000,100)
#f_pp1=real(fa_cromer('Si',Ex))
##f=elements['Si']
##f_pp1=array(f[1])+float(f[3])
##f_pp2=real(fa_asf('Si',Ex))
#sf_p=UnivariateSpline(E,array(f[1]))          # this is what we choose to interpolate
#sf_p2=KroghInterpolator(E,array(f[1]))        # these drop too rapidly at high energies
#sf_p3=BarycentricInterpolator(E,array(f[1]))  # identical to the previous one
#f_ps  = sf_p(Ex)
#f_ps2 = sf_p2(Ex)
#f_ps3 = sf_p3(Ex)
##plot(E,f_pp1, 'bo')
##plot(Ex,f_pp2, 'g-')
#plot(Ex,f_ps3, 'k-')
#plot(Ex,f_pp,'b-')
##show()

################################################################################
### STRUCTURAL SUM 
################################################################################

def expFh_dia(h,k,l):
      
      H=2.0*pi*array([h,k,l])
      
      r=[[]]*8
      r[0]=array([0.00,0.00,0.00])
      r[1]=array([0.25,0.25,0.25])
      r[2]=array([0.50,0.50,0.00])
      r[3]=array([0.75,0.75,0.25])
      r[4]=array([0.50,0.00,0.50])
      r[5]=array([0.00,0.50,0.50])
      r[6]=array([0.75,0.25,0.75])
      r[7]=array([0.25,0.75,0.75])
            
      Fh=0
      for x in r:
            Fh=Fh+exp(1j*sum(H*x))      
      return Fh
#----------------------------------------------------------
def expFh_sph(elem,h,k,l):
    
    H=2.0*pi*array([h,k,l])    
    t=[[]]*4
    t[1]=array([0.0,0.0,0.0])
    t[2]=array([2.0/3.0,1.0/3.0,1.0/3.0])
    t[3]=array([1.0/3.0,2.0/3.0,2.0/3.0])
            
    if elem=='Al':
      z=0.35220
      p=[[]]*5
      p[1]=array([0.0,0.0,z])
      p[2]=array([0.0,0.0,0.5-z])
      p[3]=array([0.0,0.0,-z])
      p[4]=array([0.0,0.0,0.5+z]) 
            
      r=[[]]*12
      for k in (1,2,3):
        for m in (1,2,3,4):
          n=m+4*(k-1)
          r[n-1]=p[m]+t[k]
      
    elif elem=='O':
      x=0.30627       
      p=[[]]*7
      p[1]=array([x,0.0,0.25])
      p[2]=array([0.0,x,0.25])
      p[3]=array([-x,-x,0.25])
      p[4]=array([-x,0.0,0.75])
      p[5]=array([0.0,-x,0.75])
      p[6]=array([x,x,0.75])

      r=[[]]*18
      for k in (1,2,3):
        for m in (1,2,3,4,5,6):
          n=m+6*(k-1)
          r[n-1]=p[m]+t[k]

    Fh=0
    for x in r:
        Fh=Fh+exp(1j*sum(H*x))      
    return Fh
# ----------------------------------------------------------------
def expFh_pyr(elem,h,k,l):   # Bailyss AM 1977
    
    H=2.0*pi*array([h,k,l])    
            
    if elem=='Fe':
      r=[[]]*4
      r[0]=array([0.00100,0.00200,0.00300])
      r[1]=array([0.49660,0.00010,0.50360])
      r[2]=array([0.50010,0.50200,0.00110])
      r[3]=array([-0.00060,0.50130,0.50380])
      
    elif elem=='S':
      r=[[]]*8
      r[0]=array([0.38570,0.38320,0.38400])
      r[1]=array([0.11490,0.61140,0.88460])
      r[2]=array([0.88540,0.11570,0.61430])
      r[3]=array([0.61530,0.88650,0.11410])
      r[4]=array([0.61510,0.61320,0.61370])
      r[5]=array([0.88540,0.38180,0.11490])
      r[6]=array([0.11470,0.88560,0.38410])
      r[7]=array([0.38570,0.11610,0.88420])
      
    Fh=0
    for x in r:
        Fh=Fh+exp(1j*sum(H*x))      
    return Fh

# ---------------------------------------------------------------------------------
def expFh_SiC4H(elem,h,k,l):   #Bauer AC 2001  +1 model agrees with the refinement
    
    H=2.0*pi*array([h,k,l])    
    delta1=0.0        ; eps1=0.0
    delta2=0.0        ; eps2=0.0
        
    if elem=='Si':
      ksi=delta1 
      tau=4.0/16.0+delta2
    elif elem=='C':
      ksi=3.0/16.00+eps1
      tau=7.0/16.0+eps2
    
    r=[[]]*4
    r[0]=array([0.0,0.0,ksi])
    r[1]=array([2.0/3.0,1.0/3.0,0.5+ksi])    
    r[2]=array([1.0/3.0,2.0/3.0,tau])  
    r[3]=array([1.0/3.0,2.0/3.0,0.5+tau])
        
    Fh=0
    for x in r:
        Fh=Fh+exp(1j*sum(H*x))      
    return Fh

# ---------------------------------------------------------------------------------
def expFh_SiC6H(elem,h,k,l):   #Bauer AC 2001  +1 model agrees with the refinement
    
    H=2.0*pi*array([h,k,l])    
    delta1=0.0 ; eps1=0.0
    delta2=0.0 ; eps2=0.0
    delta3=0.0 ; eps3=0.0
        
    if elem=='Si':
      ksi=delta1 
      tau=4.0/24.0+delta2
      v=8.0/24.0+delta3
    elif elem=='C':
      ksi=3.0/24.0+eps1
      tau=7.0/24.0+eps2
      v=11.0/24.0+eps3
    
    r=[[]]*6
    r[0]=array([0.0,0.0,ksi])
    r[1]=array([0.0,0.0,0.5+ksi])
    r[2]=array([1.0/3.0,2.0/3.0,tau])  
    r[3]=array([2.0/3.0,1.0/3.0,0.5+tau])
    r[4]=array([2.0/3.0,1.0/3.0,v])
    r[5]=array([1.0/3.0,2.0/3.0,0.5+v])
        
    Fh=0
    for x in r:
        Fh=Fh+exp(1j*sum(H*x))      
    return Fh
#------------------------------------------------------------------------------------
def expFh_SiO2(elem,h,k,l):   #Le Page et al 1980
    
    H = 2.0*pi*array([h,k,l])    
    u0 = 0.46981
    x0 = 0.41372
    y0 = 0.26769
    z0 = 0.11880  # positive because definition of atom positions from Nuttall 1981
                  # z0 < 0 ???  
        
    if elem=='Si':    # Wyckoff posiitons from www.cryst.ehu.es
      r=[[]]*3
      r[0] = array([u0,0.0,2.0/3.0])
      #r[0] = array([u0,0.0,0.0])  often in literature but not Wyckoff position that corresponds to x0 y0 z0
      r[1] = array([-u0,-u0,0.0])
      r[2] = array([0.0,u0,1.0/3.0])    
    #
    elif elem=='O':    
      r=[[]]*6
      r[0] = array([x0,y0,z0])
      r[1] = array([y0-x0,-x0,z0+1.0/3.0])
      r[2] = array([-y0,x0-y0,z0+2.0/3.0])
      #r[3] = array([x0-y0,-y0,z0])
      r[3] = array([y0,x0,-z0])
      r[4] = array([x0-y0,-y0,1.0/3.0-z0])
      r[5] = array([-x0,y0-x0,2.0/3.0-z0])
      #
      #r[4] = array([y0,x0,2.0/3.0-z0])
      #r[5] = array([-x0,y0-x0,1.0/3.0-z0])
                    
    Fh=0
    for x in r:
        Fh=Fh+exp(1j*sum(H*x))      
    return Fh


############################################################################
##  DWF  From thermal ellipsoid
############################################################################
def debye_tellipse(ind,lp,angles,U):
    from constants import r2d
    h,k,l = ind
    a,b,c = lp
    alpha,beta,gamma = angles
    alpha = alpha/r2d
    beta = beta/r2d
    gamma = gamma/r2d
    U11,U22,U33,U12,U13,U23 = U
    
    # formulas from 1.1.3 International Tables for Crystallography
    V = a*b*c*sqrt(1.0-(cos(alpha))**2.0-(cos(beta))**2.0-(cos(gamma))**2.0 + 2.0*cos(alpha)*cos(beta)*cos(gamma))
    
    a_ = b*c*sin(alpha)/V
    b_ = c*a*sin(beta)/V
    c_ = a*b*sin(gamma)/V
    
    cos_alpha_ = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma))
    cos_beta_ =  (cos(gamma)*cos(alpha)-cos(beta))/(sin(gamma)*sin(alpha))
    cos_gamma_ = (cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta))
        
    summ = U11*(h*a_)**2.0 + U22*(k*b_)**2.0 + U33*(l*c_)**2.0 + \
           2.0*U12*h*k*a_*b_*cos_gamma_ + \
           2.0*U13*h*l*a_*c_*cos_beta_  + \
           2.0*U23*k*l*b_*c_*cos_alpha_ 
    
    return exp(-2.0*pi**2.0*summ)

############################################################################
##  DWF  Sears and Sheley Acta Cryst. A47,441 (1991)
############################################################################
def debye_sears(atom,T,qx):
    
    if atom=='C':
      M=12.011
      nu_m=40.09
      #Tm=1924.0
      Tm=2200.0
      alp=2.352
      f_2=3.114
      f_1=1.552
      f2=0.576
    elif atom=='Si':
      M=28.086
      nu_m=15.85
      Tm=761.0
      alp=6.933
      f_2=6.943
      f_1=2.136
      f2=0.497
    elif atom=='Ge':
      M=72.59
      nu_m=9.20
      Tm=442.0
      alp=4.518       
      f_2=7.495       
      f_1=2.245       
      f2=0.476       
#----------------------------
    elif atom=='Al':
      M=26.982
      nu_m=9.75
      Tm=468
      alp=3.670
      f_2=4.029
      f_1=1.778
      f2=0.445
#----------------------------
    elif atom=='Fe':
      M=55.847
      nu_m=9.54
      Tm=458
      alp=3.079
      f_2=3.310
      f_1=1.602
      f2=0.522
#----------------------------
    elif atom=='Cu':
      M=63.546
      nu_m=7.29
      Tm=350.0
      alp=3.681
      f_2=3.737
      f_1=1.699
      f2=0.479
                    
    y=T/Tm  #t=T/Tm
#    Jy=[]
#    for y in t:
    if y<0.2:    
        Jy=f_1+pi**2.0/3.0*alp*y**2.0
    else:
        Jy=2.0*f_2*y+1.0/(6.0*y)-f2/(360.0*y**3.0)
        
    B=39.904/(M*nu_m)*Jy
    W=B*qx**2.0
    return [B,exp(-W)]
