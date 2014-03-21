#!/usr/bin/env python
########################################################
## SET OF FUNCTIONS TO CALCULATE STRUCTURE FACTORS
########################################################
# VER=0.17
#########################################################################
# to do list:
# 1. check alternatives for DWF
#########################################################################

import dtxrd
import os
libpath = os.path.dirname(dtxrd.__file__)
from numpy import *
#
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
###########################################################################
###  ANOMALOUS CORRECTIONS
###########################################################################
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

def fa_asf(atom,Ex):
    #import xdp
    from myio import readFile
    fileName=libpath+'/asf/'+atom+'.asf'
    d1,d2=readFile(fileName)
    Easf=d2[:,0]; Easf=1.0e3*Easf
    f1=d2[:,1]
    f2=d2[:,2]
#    data=open(fileName, 'r')
#    line0=data.readline(); print line0
#    while 1:
#        line=data.readline()
#        stuff=line.split(' ')
#        if not(line):
#          break
    Z=float(elementsZ[atom])
    f_1 = interp1d(Easf,f1)  # it seems that f1 in fa_asf is f1=f_p+f_rel+Z
    f_2= interp1d(Easf,f2)

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
    return exp(-W)
