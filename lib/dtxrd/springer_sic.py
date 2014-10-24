#Thermal expansion of SiC-4H and SiC-6H
#Springer http://link.springer.com/chapter/10.1007/10832182_600#page-1
#######################################################################
from numpy import *

# original data for SiC-4H from Z Li and R.C. Bradt J Appl. Phys. 60, 612 (1986)
def a_SiC4H(T):        
    return 3.0801 + 9.9031e-6*T + 5.4833e-9*T**2.0 - 1.6613e-12*T**3.0 # Angstrom
    
def c_SiC4H(T):
    return 10.085 + 3.1124e-5*T + 1.3266e-8*T**2.0 - 3.6252e-12*T**3.0 # Angstrom

# presumably from Taylor and Jones "Silicon Carbide" book
def a_SiC6H(T):
    return 10.0*(0.30813 + 1.0064e-6*T + 5.0126e-10*T**2.0 - 1.3982e-13*T**3.0) # Angstrom
    
def c_SiC6H(T):
    return 10.0*(1.5116 +4.8006e-6*T + 1.8754e-10*T**2.0 - 4.2862e-13*T**3.0) # Angstrom
            

        
    
    
    
