#!/usr/bin/python

from curvestat import lorentz

#from numpy import *
#from pylab import *

def CuKa(e,w,q,x):
    
  result=0
  for (ex,wx,qx) in map(None,e,w,q):
        b=[0.0,qx,ex,0.5*wx]
        result=result+lorentz(b,x)
  return result

# Hartwig JAC1993 single crystal
def CuKa_hartwig_sc(x):
  e=[8047.838,8045.401,8028.022,8026.467]
  w=[2.26,3.38,2.71,3.44]
  q=[0.956,0.092,0.330,0.113]
  return CuKa(e,w,q,x)

# Hartwig JAC1993 double crystal
def CuKa_hartwig_dc(x):
  e=[8047.839,8045.343,8027.939,8026.351]
  w=[2.24,3.71,2.86,3.25]
  q=[0.950,0.098,0.382,0.082]
  return CuKa(e,w,q,x)

# Berger XS1986
def CuKa_berger_ss(x):
  e=[8047.838,8045.425,8027.963,8026.347]
  w=[2.28,3.1,2.70,3.2]
  q=[0.970,0.095,0.331,0.093]
  return CuKa(e,w,q,x)


def MoKa(x):
  e=[17479.3,17374.3]
  w=[6.31,6.49]
  q=[1.0,0.5]
  return CuKa(e,w,q,x)

#TEST
#ev=arange(8020,8060,0.1)
#har_sc=CuKa_hartwig_sc(ev)
#har_dc=CuKa_hartwig_dc(ev)
#ber_ss=CuKa_berger_ss(ev)
#figure()
#plot(ev,har_sc,'b-')
#plot(ev,har_dc,'r-')
#plot(ev,ber_ss,'g-')
#show()










        


    