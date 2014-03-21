#!/usr/bin/env python

'''
a subroutine to calculate parameters of a Bragg reflection

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

################################################################
# v 0.02 #######################################################
################################################################

from numpy import *
from okada_si import *
from stoupin_c import *
from carr_ge import *
from lucht_sph import *
from fh import *

from dtxrd0 import *
from dtxrd2_k import *
from constants import *
from curvestat import *
from chi import *
#----------------------------------------------------------------------
# Constants
#----------------------------------------------------------------------
#r2d=180.0/pi        # radians to degrees conversion
#----------------------------------------------------------------------
#hpl=4.13566733e-15     # [eV*s] Planck constant
#hpl_bar=1.05457266e-34 # [J*s] Planck constant
#qe=1.60217653e-19     # [C] Electron charge
#me=9.0193897e-31     # [kg] mass of electron
#cl=299792458.0e10      # [Angstrom/s] speed of light
#re=2.8179402894e-5  #[A] classical radius of electron
#----------------------------------------------------------------------

def thc_find(Ex,eta,phi,dc,crystal,P):
    dh=crystal[1]; Eb=0.5*hpl*cl/dh
    thb=arcsin(Eb/Ex)
    thc=thb; ee=1.0; i=0; imax=100
    #
    while (ee > 1.0e-6):
        i=i+1                       #; print i
        if i > imax:
            print ('thc_find convergence problem!')
            break
        else:
            k=[cos(thc)*cos(phi),cos(thc)*sin(phi),-sin(thc)]
            k_pr,thx,dth,eps,bh,Ty,Ry=dtxrd0_k(k,eta,dc,Ex,P,crystal)
            ee=abs(thx-thc)/dth      #; print eps
            thc=thx
            if 1.0-P < 1.0e-6:
                P=cos(2.0*thx)
            thc_pr=arcsin(k_pr[2])

    return [thc,thc_pr,eps,bh,Ty,Ry]

################################################################################
## FIND TH_MAX
################################################################################

def thmax_find(Ex,eta,phi,dc,crystal,P,flag):

    thc,thc_pr,eps,bh,Ty,Ry=thc_find(Ex,eta,phi,dc,crystal,P)
#  Chi,wh_s,wh,eps_s,eps,dth_s,dth,de,Tplot,Rplot,dth_pr,eps_pr
    result0=dtxrd1(thc,eta,phi,dc,Ex,P,crystal)
#   Chi=result0[0]
    [dth_s,dth,dth_pr]=result0[3]
    thv=thc+1.5*dth*arange(-1.0,1.0,2.0/1000.0)

    T=[]; R=[]
    for thx in thv:
        k=[cos(thx)*cos(phi),cos(thx)*sin(phi),-sin(thx)]
        [k_pr,Tx,Rx]=dtxrd1_k(k,eta,dc,Ex,P,crystal)
        T=T+[Tx]; R=R+[Rx]

    thv=array(thv); R=array(R); T=array(T)
    if flag=='R':   stat1=curvestat(thv,R,0.0)
    elif flag=='T': stat1=curvestat(thv,T,0.0)
    else: stat1=0; print " flag input error "

    thmax=stat1[0]
    k=[cos(thmax)*cos(phi),cos(thmax)*sin(phi),-sin(thmax)]
    [k_pr,Tx,Rx]=dtxrd1_k(k,eta,dc,Ex,P,crystal)
    thmax_pr=arcsin(k_pr[2])

    return [thmax,thmax_pr,eps,bh,Ty,Ry]

#################################################################################
## FIND TH_R
#################################################################################
def thr_find(el,h,k,l,eta,phi,T,dc,Ex):
#def thr_find(Ex,eta,phi,dc,crystal,P):
    eta=1.0e-6+eta
    E0=Ex; ee=1.0; i=0; imax=100
    while (ee > 1.0e-16):
        i=i+1
        if i > imax:
            print ('thr_find convergence problem!')
            break
        else:
            [[Chi0,Chih,Chih_],dh],flagFh=chi(el,h,k,l,T,E0)
            #
            lam0=hpl*cl/E0
            wh_s=-2.0*real(Chi0)*(dh/lam0)**2.0  #; print 'wh_s = ', wh_s
            # theta_r:
            denr=1.0+sqrt(1.0+4.0*wh_s*(tan(eta))**2.0)
            bthr=abs(arctan(2.0*wh_s*tan(eta)/denr))
            thr=0.5*pi-bthr
            # lambda_r:
            lamr=2.0*dh*cos(bthr)*(1.0-tan(bthr)/tan(eta))
            Er=hpl*cl/lamr
            ee=abs(E0-Er)/Er
            E0=Er
#   print " thr_find # of iterations = ", i
#   print "wh_s = ", wh_s
    return [thr,Er]
#################################################################################
