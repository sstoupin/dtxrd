#Thermal expansion of alpha-quartz
#Barron JPC1982
####################################
from numpy import *
from scipy import integrate
from scipy.integrate import simps
#from scipy.optimize import *
from scipy.interpolate import interp1d


# These are only for low temperatures
#def alpha_a(x):
#  return x**3.0*(12.4 + 5.5*(x/10.0)**2.0 + 6.7*(x/10.0)**4.0 + 2.2*(x/10.0)**6.0)*1.0e-12 #K-1
  
#def alpha_c(x):
#  return (-14.5*x**3.0 + 0.073*x**5.0)*1.0e-12 #K-1

t_dat = array([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,35,40,50,57.5,65,75,85,100,125,150,175,200,225,250,273,283,293,298,373,473,573,673,773,813,833,838,843])
a_dat = array([0.18,0.33,0.59,0.99,1.65,2.70,4.26,6.47,9.48,13.29,17.95,23.48,29.85,37.09,45.11,53.7,72.5,93.0,114.7,137.2,160,217,272,373,441,500,566,631,720,846,960,1051,1134,1207,1272,1329,1352,1375,1387,1555,1800,2080,2580,3960,5500,7900,9500,14000])
c_dat = array([-0.16,-0.25,-0.39,-0.50,-0.60,-0.61,-0.40,-0.05,0.53,1.47,2.76,4.39,6.32,8.61,11.24,14.2,20.7,27.7,35.0,42.8,50.9,72,92,129,157,183,217,252,302,382,454,517,574,628,677,720,737,754,762,872,1036,1238,1520,2250,3250,4700,5600,8400])

t_dat = array(t_dat, dtype=float)
a_dat = array(a_dat, dtype=float)
c_dat = array(c_dat, dtype=float)

def alpha_a(x):
     alpha_a_f = interp1d(t_dat,a_dat)
     return 1.0e-8*alpha_a_f(x) 
def alpha_c(x):
     alpha_c_f = interp1d(t_dat,c_dat)
     return 1.0e-8*alpha_c_f(x) 

a0 = 4.9130       # Cohen and Sumner /Keith 1950
c0 = 5.4046
T0 = 298.0

def a_sio2(x):
  y=integrate.quad(lambda t: alpha_a(t),T0,x)
  a=a0*(1.0+y[0])
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return a

def c_sio2(x):
  y=integrate.quad(lambda t: alpha_c(t),T0,x)
  c=c0*(1.0+y[0])
#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
  return c
