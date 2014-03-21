#Thermal expansion of Ge
#Carr et al Phil. Mag. 1965
####################################
from numpy import *
from scipy import integrate
from scipy.optimize import *
from scipy.interpolate import interp1d

def alpha_carr(x):
    # 8-130 - Carr (1965)
    # 150-280 - Gibbons (1958)
    # 293 and above Touloukian v12
    data = [
    [8.0,0.35],
    [10.0,0.45],
    [12.0,0.60],
    [14.0,0.25],
    [16.0,-0.45],
    [18.0,-1.35],
    [20.0,-2.3],
    [22.0,-3.6],
    [24.0,-4.9],
    [26.0,-5.6],
    [28.0,-6.1],
    [30.0,-7.1],
    [35.0,-6.2],
    [40.0,-2.4],
    [50.0,16.5],
    [60.0,52.3],
    [70.0,91.8],
    [80.0,138.0],
    [90.0,185.0],
    [100.0,228.0],
    [110.0,269.0],
    [120.0,307.0],
    [130.0,341.0],
    [150.0,412.0],
    [200.0,482.0],
    [250.0,532.0],
    [280.0,559.0],
    [293.0,570.0],
    [400.0,620.0],
    [500.0,650.0],
    [600.0,670.0],
    [700.0,700.0],
    [800.0,720.0],
    [900.0,740.0],
    [1000.0,760.0],
    [1200.0,800.0]]

    data=transpose(array(data))
    t0=data[0]
    alp0=1.0e-8*data[1]
    alpha=interp1d(t0,alp0)
#  print " Temperature = ", x, " K"
#  print " Thermal expansion coeff. = ", alpha, " K-1"
    result=float(alpha(x))
    return result

#a0=5.6575   # Singh AC (1968)
#T0=293.15
#a0=5.6579060 # Baker AC (1975)
a0=5.657820   # Hom  AC  (1975)
T0=298.15

def a_carr(x):
    y=integrate.quad(lambda t: alpha_carr(t),T0,x)
    a=a0*(1.0+y[0])

#  print " integration abs. error", y[1]
#  print " Temperature = ", x, " K"
#  print " Lattice parameter = ", a, " Angstrom"
    return a
