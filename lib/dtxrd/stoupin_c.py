#Thermal expansion of C
#Stoupin and Shvyd'ko (2010)
####################################
from numpy import *
from scipy import integrate
from scipy.optimize import *

b=3.6e-14
c=1.21e-11
x0=212.0
dx=47.0


def alpha_stoupin(x):
    e=exp((x-x0)/dx)
    w=1.0/(1+e)
    alpha=abs(b)*x**3.0*w+abs(c)*x**2.0*(1-w)

#  print " Temperature = ", x, " K"
#  print " Thermal expansion coeff. = ", alpha, " K-1"
    return alpha

a0=3.56712
T0=298.0

def a_stoupin(x):

    y=integrate.quad(lambda t: alpha_stoupin(t),T0,x)
    a=a0*(1.0+y[0])

#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
    return a
