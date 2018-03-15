#Thermophysical properties of SiC polymorphs
from numpy import *
from scipy import integrate
from scipy.optimize import *
from scipy.interpolate import interp1d

########################################################
## Reeber & Wang (2000) MRS
########################################################

def eins(tt):
    return tt**2.0*exp(tt)/(exp(tt)-1.0)**2.0
    
def alpha_a_SiC6H_reeber(T):
    th = array([150.0,700.0,1952.5])               # K 
    x = array([2.726,40.116,6.963])*1.0e-7 # K-1
    tt = th/T
    return sum(x*eins(tt))

def alpha_c_SiC6H_reeber(T):
    th = array([125.0,600.0,1412.5])               # K 
    x = array([-12.725,59.148,3.532])*1.0e-7 # K-1
    tt = th/T
    return sum(x*eins(tt))
    
#a0 = 3.08084
#c0 = 15.11765
#T0 = 298.0

def a_SiC6H_reeber(x):
    a0 = 3.08084
    T0 = 298.0
    y=integrate.quad(lambda t: alpha_a_SiC6H_reeber(t),T0,x)
    a=a0*(1.0+y[0])
    return a

def c_SiC6H_reeber(x):
    c0 = 15.11765
    T0 = 298.0
    y=integrate.quad(lambda t: alpha_c_SiC6H_reeber(t),T0,x)
    c=c0*(1.0+y[0])
    return c

#######################################################################################
#######################################################################################

def alpha11_SiC4H_RS(T):                 # Stockmeier et al using Reeber's formula
    th11 = array([441.9,659.1,1970.01])               # K 
    x11 = array([35.5,-4.6,28.6])*1.0e-7 # K-1
    tt11 = th11/T    
    return sum(x11*eins(tt11))
def alpha33_SiC4H_RS(T):
    th33 = array([436.9,660.6,1975.1]) #K
    x33 = array([44.5,-16.7,31.1])*1.0e-7 #K-1
    tt33 = th33/T
    return sum(x33*eins(tt33))
    
def a_SiC4H_RS(x):
    #a0 = 3.07976  #RS
    a0 = 3.08051  # Bauer ACA2000 
    T0 = 300.0       
    y=integrate.quad(lambda t: alpha11_SiC4H_RS(t),T0,x)
    a=a0*(1.0+y[0])
    return a
def c_SiC4H_RS(x):
    #c0 = 10.08196 #RS
    c0 = 10.08480 # Bauer ACA2000
    T0 = 300.0
    y=integrate.quad(lambda t: alpha33_SiC4H_RS(t),T0,x)
    c=c0*(1.0+y[0])
    return c

#-------- SiC-6H -----------------------------------------------------------------------
# Stockmeier et al using Reeber's formula - reported values for the undoped sample
def alpha11_SiC6H_RS(T):                 
    th11 = array([442.4,656.4,1963.0])   # K 
    x11 = array([39.0,-5.8,29.0])*1.0e-7 # K-1
    tt11 = th11/T    
    return sum(x11*eins(tt11))

def alpha33_SiC6H_RS(T):
    th33 = array([439.0,657.0,1965.1]) #K
    x33 = array([47.9,-14.2,29.2])*1.0e-7 #K-1
    tt33 = th33/T
    return sum(x33*eins(tt33))
    
def a_SiC6H_RS(x):
    a0 = 3.08049 
    T0 = 300.0       
    y=integrate.quad(lambda t: alpha11_SiC4H_RS(t),T0,x)
    a=a0*(1.0+y[0])
    return a
def c_SiC6H_RS(x):
    c0 = 15.11508
    T0 = 300.0
    y=integrate.quad(lambda t: alpha33_SiC4H_RS(t),T0,x)
    c=c0*(1.0+y[0])
    return c

###########################################################################################
# Nakabayashi et al. MRS 2006 - low temp data
###########################################################################################
def alpha11_SiC4H_N(T):
    return (-2.0404 + 1.9374e-2*T - 1.1385e-5*T**2.0)*1.0e-6
def alpha33_SiC4H_N(T):
    return (-1.9755 + 1.8967e-2*T - 1.0971e-5*T**2.0)*1.0e-6
    

