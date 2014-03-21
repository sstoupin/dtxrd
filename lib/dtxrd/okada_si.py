#Thermal expansion of Si
#Okada and Tokumaru (1984)
####################################
from numpy import *

T1=124.0

C1=3.725e-6
C2=5.88e-3
C3=5.548e-10



def alpha_okada(T):


    alpha=C1*(1.0-exp(-C2*(T-T1)))+C3*T

    #print " Temperature = ", T, " K"
    #print " Thermal expansion coeff. = ", alpha, " K-1"

    return alpha


a0=5.430741
T0=273.2

def a_okada(T):

    a=a0+a0*(C1*(T-T0+(exp(-C2*(T-T1))-exp(-C2*(T0-T1)))/C2)+0.5*C3*(T**2-T0**2))

    #print " Temperature = ", T, " K"
    #print " Lattice parameter = ", a, " Angstrom"

    return a
