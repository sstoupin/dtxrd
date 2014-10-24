#Thermal expansion of pyrite 
#R.S. B. Chrystall 1965
####################################
from numpy import *
from scipy import integrate
from scipy.optimize import *

def alpha_chrystall(x):
  x = x - 273.15 # convert Kelvin to Celcius
  alpha = 8.87e-6 + 5.82e-9*x
#  print " Temperature = ", x, " K"
#  print " Thermal expansion coeff. = ", alpha, " K-1"  
  return alpha

a0=5.4174       # Chrystall 1965
T0=273.15+16.0

def a_chrystall(x):
  y=integrate.quad(lambda t: alpha_chrystall(t),T0,x)
  a=a0*(1.0+y[0])
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return a

