#Thermal expansion of Sapphire alpha-Al2O3
#Lucht et al (2003), Shvyd'ko et al (2002)
############################################
from numpy import *

def x_lucht(x04,x44,x01,x11,tx,dtx,T):

    w=(1.0+exp((sqrt(T)-sqrt(tx))/sqrt(dtx)))**(-1.0)
    x=(x44*T**4.0+x04)*w+(x11*T+x01)*(1.0-w)
    return x

def a_lucht(T):

    a04=4.756274  #Angstrom
    a44=1.5e-12   #Angstrom*K-4
    a01=4.7501    #Angstrom
    a11=2.95e-5   #Angstom*K-1
    ta=193        #K
    dta=1.0       #K
    result=x_lucht(a04,a44,a01,a11,ta,dta,T)
    return result

def c_lucht(T):

    c04=12.981943  #Angstrom
    c44=5.8e-12    #Angstrom*K-4
    c01=12.9633    #Angstrom
    c11=9.2e-5     #Angstom*K-1
    tc=187         #K
    dtc=0.99       #K
    result=x_lucht(c04,c44,c01,c11,tc,dtc,T)
    return result
