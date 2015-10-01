###########################################################
## Calculate numerical derivative
# difference method + interpolation
# a - xvector
# b - yvector
# n - spline points
###########################################################
from scipy import array, linspace
from scipy.interpolate import UnivariateSpline, interp1d

def derivs(a,b,n):
    x=0.5*(a[:len(a)-1]+a[1:])
    dx=0.5*(a[1:]-a[:len(a)-1])
    dy=0.5*(b[1:]-b[:len(a)-1])        
    dery=dy/dx
    xmin=x[0]; xmax=x[len(x)-1]
    xs=linspace(xmin,xmax,n)
    s=UnivariateSpline(x,dery,s=1)
    ys=s(xs)
    return [xs,ys]
    
def deriv(a,b):
    x=0.5*(a[:len(a)-1]+a[1:])
    dx=0.5*(a[1:]-a[:len(a)-1])
    dy=0.5*(b[1:]-b[:len(a)-1])        
    dery=dy/dx
    return [x,dery]        

def derivi(x0,y0,n):
    f_int=interp1d(x0,y0) 
    xmin=x0[0]; xmax=x0[len(x0)-1]
    a = linspace(xmin,xmax,n)
    b = f_int(a)
    x=0.5*(a[:len(a)-1]+a[1:])
    dx=0.5*(a[1:]-a[:len(a)-1])
    dy=0.5*(b[1:]-b[:len(a)-1])        
    dery=dy/dx
    return [x,dery]        
    