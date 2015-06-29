#Thermophysical properties of Si
####################################
from numpy import *
from scipy import integrate
from scipy.optimize import *
from scipy.interpolate import interp1d

# Okada & Tokumaru JAP 1984 
# Good only from 124 to 1500 K (Shvyd'ko book)
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


### REEBER  Phys. Stat. Sol (a) 1975
###
def eins(tt):
    return tt**2.0*exp(tt)/(exp(tt)-1.0)**2.0
def alpha_reeber(T):
    th = array([30.0, 68.75, 189.4, 664.4])               # K 
    x = array([-0.02953, 0.23097, -8.707, 50.949])*1.0e-7 # K-1
    tt = th/T
    return sum(x*eins(tt))
        
def alpha2_reeber(T):
    result = []
    th = array([30.0, 68.75, 189.4, 664.4])               # K 
    x = array([-0.02953, 0.23097, -8.707, 50.949])*1.0e-7 # K-1
    for t in T:
        tt = th/t
        result = result + [sum(x*eins(tt))]
    return array(result)

def a_reeber(x):
    y=integrate.quad(lambda t: alpha_reeber(t),T0,x)
    a=a0*(1.0+y[0])
    #  print " integration abs. error", y[1]
    #  print " Temperature = ", x, " K"
    #  print " Lattice parameter = ", a, " Angstrom"
    return a
    
def a2_reeber(T):
    result = []
    for x in T:    
       y=integrate.quad(lambda t: alpha_reeber(t),T0,x)
       a=a0*(1.0+y[0])
       result = result + [a]
    #  print " integration abs. error", y[1]
    #  print " Temperature = ", x, " K"
    #  print " Lattice parameter = ", a, " Angstrom"
    return array(result)
#                
################################################################################################################
# Thermal conductivity
# Glassbrenner and Slack PR 1964
# T[K], k_Si [W cm-1 K-1]

data = [
        [4.381866938,  2.12340907972],
        [5.18722025277,  4.12956957652],
        [7.45569860518,  8.04352190764],
        [10.1803716825,  16.155634639],
        [14.2246812416,  22.8101537211],
        [15.418018613,   26.6669935896],
        [17.7648743112,  31.8449130737],
        [19.2611123992,  35.7165049011],
        [25.6417315045,  35.0629416865],
        [32.4415407668,  33.7010055799],
        [43.6652707228,  29.8273822809],
        [66.634307202,   17.6337813104],
        [92.7694514846,  10.2035136567],
        [138.918905786,  4.90145144886],
        [195.510711425,  2.61054010128],
        [301.463390863,  1.51177612031],
        [410.932592999,  0.960194855556],
        [484.515617326,  0.806031103915],
        [559.637728842,  0.690688371369],
        [626.757026128,  0.610406798998],
        [653.844837035,  0.501394925644],
        [709.186451581,  0.533934055688],
        [717.016471088,  0.481371647611],
        [739.496481117,  0.466738506262],
        [762.739771257,  0.44788087418],
        [811.441010287,  0.412420541054],
        [845.73212756,  0.383664605058],
        [899.732461687,  0.353288503982],
        [909.736025045,  0.315223183704],
        [947.308985197,  0.332108869421],
        [1007.79504225,  0.305814620609],
        [1050.94793245,  0.264569645141],
        [1105.33329742,  0.287574476657],
        [1129.00963263,  0.25928566351],
        [1151.86750499,  0.273130541016],
        [1201.00525394,  0.241246218409],
        [1212.59058395,  0.262137837673],
        [1331.07196954,  0.219925288438],
        [1344.22115251,  0.231649650678],
        [1428.2968231,  0.251812482284],
        [1536.03506956,  0.213430454321],
        [1551.09006992,  0.227152283951]]
data = transpose(data)
temp = data[0]
k0 =data[1]
                
def k_si(x):
    k_=interp1d(temp,k0)    
    return k_(x)  # units W cm-1 K-1
        
    
    