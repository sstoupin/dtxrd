#Thermal expansion of C
####################################
from numpy import *
from scipy import integrate
from scipy.optimize import *

# Stoupin & Shvyd'ko PRB 2011
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
  #a1=a0*(1.0+y[0])
  a2=a0*exp(y[0])
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return a2

# vector format
def av_stoupin(temp):
  stuff = []
  for x in temp:
    y=integrate.quad(lambda t: alpha_stoupin(t),T0,x)
    #a1=a0*(1.0+y[0])
    a2=a0*exp(y[0])
    stuff = stuff +[a2]
  stuff = array(stuff)  
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return stuff
  
###### New Stuff #################################
##################################################

def eins(tt):
    return tt**2.0*exp(tt)/(exp(tt)-1.0)**2.0

def alpha_reeber(T):
    th = array([200.0,880.0,2137.5])
    x = array([0.4369e-7,15.7867e-7,42.5598e-7])
    tt = th/T
    return sum(x*eins(tt))

def alpha_jacobson(T):
    th = array([302.8,1007.5,2249.0])
    x = array([0.6339e-7,20.2700e-7,38.0149e-7])
    tt = th/T
    return sum(x*eins(tt))

def alpha_WK(T): # Widl & Koidl 2001 textbook
    T = T - 273.15
    return 8.19e-7 + 1.107e-8*T - 1.48e-11*T**2.0 + 1.08e-14*T**3.0
    
def alpha_JRWS1(T):
    x = 1.0e-6*array([0.0118, 0.2932, 2.9149, 2.1331]) 
    th = array([176.8, 568.8, 1273.3, 2250.5])
    tt = th/T
    return sum(x*eins(tt))

def alpha_JRWS2(T):
    x = 1.0e-6*array([0.025, 0.4212, 3.4969, 2.3454]) 
    th = array([240.6, 655.2, 1382.6, 3232.8])
    tt = th/T
    return sum(x*eins(tt))

def av_JRWS1(temp):
  stuff = []
  for x in temp:
      y=integrate.quad(lambda t: alpha_JRWS1(t),T0,x)
      #a1=a0*(1.0+y[0])
      a2=a0*exp(y[0])
      stuff = stuff +[a2]
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return array(stuff)

def av_JRWS2(temp):
  stuff = []
  for x in temp:
      y=integrate.quad(lambda t: alpha_JRWS2(t),T0,x)
      #a1=a0*(1.0+y[0])
      a2=a0*exp(y[0])
      stuff = stuff + [a2]
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return array(stuff)

#################################################
def alpha_JS1(T):
    x = 1.0e-6*array([0.0096, 0.2656, 2.6799, 2.3303])    
    th = array([159.3,548.5,1237.9,2117.8])
    tt = th/T
    return sum(x*eins(tt))

def alpha_JS2(T):
    x = 1.0e-6*array([0.021000, 0.389700, 3.444700, 2.279600]) 
    th = array([225.2, 634.0, 1364.5, 3068.8])
    tt = th/T
    return sum(x*eins(tt))

def a_JS1(x):
    y=integrate.quad(lambda t: alpha_JS1(t),T0,x)
    #a1=a0*(1.0+y[0])
    a2=a0*exp(y[0])
    return a2

def a_JS2(x):
    y=integrate.quad(lambda t: alpha_JS2(t),T0,x)
    #a1=a0*(1.0+y[0])
    a2=a0*exp(y[0])
    return a2

def av_JS1(temp):
  stuff = []
  for x in temp:
      y=integrate.quad(lambda t: alpha_JS1(t),T0,x)
      #a1=a0*(1.0+y[0])
      a2=a0*exp(y[0])
      stuff = stuff +[a2]
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return array(stuff)

def av_JS2(temp):
  stuff = []
  for x in temp:
      y=integrate.quad(lambda t: alpha_JS2(t),T0,x)
      #a1=a0*(1.0+y[0])
      a2=a0*exp(y[0])
      stuff = stuff + [a2]
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return array(stuff)

