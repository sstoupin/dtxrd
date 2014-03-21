#!/usr/bin/python
# Wrapper function for leastsq (scipy)
# ------------------------------------------------------------ 
# Author: Stanislav Stoupin (sstoupin@aps.anl.gov) 2011
#-------------------------------------------------------------
# version 0.2 

from numpy import *
from scipy.optimize import *

def fit1d(residuals, p_init, x, y):
      
      N=len(x)  
      fit=leastsq(residuals, p_init, args=(y,x), full_output=1, ftol=1e-32, xtol=1e-32)
      # returns:
      # [0] The solution (or the result of the last iteration for an unsuccessful call).
      # [1] Uses the fjac and ipvt optional outputs to construct an estimate of the jacobian 
      # around the solution. None if a singular matrix encountered (indicates very flat 
      # curvature in some direction). This matrix must be multiplied by the residual standard 
      # deviation to get the covariance of the parameter estimates - see curve_fit.
      p=fit[0]
      dte=residuals(p,y,x)
            
      chisq = abs(sum(dte**2.0))   
      r2 = sum(dte**2.0)/sum(y**2.0)   # look into IXS error standards
      sstd = sqrt(1.0/N*chisq)             # this is just definition of std
      stat=[chisq,r2,sstd]        

      if fit[1] == None:
          print fit[3]
          l1=[]
          for x in p:
              l1.append('N/A')
              uncert=l1                       
          
          return [p,uncert,dte,stat,'NA']
                   
      else:    
          corr = sstd*sqrt(abs(fit[1])) # multiply by std to get the covariance matrix          
          # square root is taken from the cov_x to convert into 
          # author's comment: the conv_x appears to be normalized on dispersion for the data (std^2) 
          uncert = diag(corr)

          return [p,uncert,dte,stat,sstd**2.0*fit[1]]
      

      