#Thermal expansion of SiC-4H and SiC-6H
#Springer http://link.springer.com/chapter/10.1007/10832182_600#page-1
#######################################################################
from numpy import *
from scipy import integrate
from scipy.optimize import *
from scipy.interpolate import interp1d


# original data for SiC-4H from Z Li and R.C. Bradt J Appl. Phys. 60, 612 (1986)
def a_SiC4H(T):        
    T = T - 273.15 # convert Kelvin to Celcius
    return 3.0801 + 9.9031e-6*T + 5.4833e-9*T**2.0 - 1.6613e-12*T**3.0 # Angstrom
###
def alpha_a_SiC4H(T):        
    T = T - 273.15 # convert Kelvin to Celcius
    #return 3.0801 + 9.9031e-6*T + 5.4833e-9*T**2.0 - 1.6613e-12*T**3.0 # Angstrom
    return 3.21e-6 + 3.56e-9*T - 1.62e-12*T**2.0
### 

T0 = 298.0
a0 = 3.073
c0 = 10.053

def a_SiC4H_int(x):
    y=integrate.quad(lambda t: alpha_a_SiC4H(t),T0,x)
    a=a0*(1.0+y[0])
    #  print " integration abs. error", y[1]
    #  print " Temperature = ", x, " K"
    #  print " Lattice parameter = ", a, " Angstrom"
    return a
   
def c_SiC4H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    return 10.085 + 3.1124e-5*T + 1.3266e-8*T**2.0 - 3.6252e-12*T**3.0 # Angstrom
def alpha_c_SiC4H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    #return 10.085 + 3.1124e-5*T + 1.3266e-8*T**2.0 - 3.6252e-12*T**3.0 # Angstrom
    return 3.09e-6 + 2.63e-9*T - 1.08e-12*T**2.0
    
def c_SiC4H_int(x):
    y=integrate.quad(lambda t: alpha_c_SiC4H(t),T0,x)
    c=c0*(1.0+y[0])
    #  print " integration abs. error", y[1]
    #  print " Temperature = ", x, " K"
    #  print " Lattice parameter = ", a, " Angstrom"
    return c
    
###        
# presumably from Taylor and Jones "Silicon Carbide" book
# from Li&Bradt JACerSoc 1987
def a_SiC6H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    return 10.0*(0.30813 + 1.0064e-6*T + 5.0126e-10*T**2.0 - 1.3982e-13*T**3.0) # Angstrom
def alpha_a_SiC6H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    #return 10.0*(0.30813 + 1.0064e-6*T + 5.0126e-10*T**2.0 - 1.3982e-13*T**3.0) # Angstrom
    return 3.27e-6 +3.25e-9*T - 1.36e-12*T**2.0
        
def c_SiC6H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    return 10.0*(1.5116 +4.8006e-6*T + 1.8754e-10*T**2.0 - 4.2862e-13*T**3.0) # Angstrom
def alpha_c_SiC6H(T):
    T = T - 273.15 # convert Kelvin to Celcius
    #return 10.0*(1.5116 +4.8006e-6*T + 1.8754e-10*T**2.0 - 4.2862e-13*T**3.0) # Angstrom
    return 3.18e-6 + 2.48e-9*T - 8.51e-13*T**2.0
    
            

