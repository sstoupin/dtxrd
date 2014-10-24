#!/usr/bin/env python

'''
a subroutine to calculate yield of secondary radiation in dynamical diffraction condition

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

from numpy import *
from deriv import *
from dtxrd0 import *
#----------------------------------------------------------------------------------------------
# v 0.01
#----------------------------------------------------------------------------------------------

def syield([GG,DD],zv,Ls):        
    [GG,DD]=stuff
    G0=GG[0]; Gh=GG[1]
    D0=DD[0]; Dh=DD[1]
    D02=abs(D0)**2.0; [zz,dD02]=deriv(zv,D02)
    Dh2=abs(Dh)**2.0; [zz,dDh2]=deriv(zv,Dh2)
    kappaz=-G0*dD02-Gh*dDh2   # -Gh gets positive when Gh<0 (Bragg) and negative when Gh>0 (Laue)
    Pz=exp(-zz/Ls)
    syield_tot=zstep*sum(syieldz)
     
   
    
    

