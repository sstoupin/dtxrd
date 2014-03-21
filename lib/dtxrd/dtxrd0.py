#!/usr/bin/env python

'''
a subroutine to calculate a set of parameters for transmitted/reflected a monochromatic wave

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

#from numpy import *
from fh import *
from constants import *
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

def dtxrd(element,h,k,l,thb,eta,phi,a,dh,T,dc,Ex,P):

###############################################################################
##  GEOMETRY
###############################################################################
    lamx=hpl*cl/Ex
    K0=2.0*pi/lamx   # ; print "K0 = ", K0
    H0=2.0*pi/dh     # ; print "H0 = ", H0
    G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
    Gh=G0-H0/K0*cos(eta)                            #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                    #; print "Gh = ", Gh # same actually
    bh=G0/Gh
###############################################################################
## STRUCTURE FACTOR and SUSCEPTIBILITIES
###############################################################################
########
    qx=0.5/dh #; print "qx (dh) = ", qx
#        qx1=sin(thb)/(hpl*cl/Ex); print "qx (lam) = ", qx1
    f0h=f0h_ICD(element,qx) #; print "f0h = ", f0h
    f00=f0h_ICD(element,0)  #; print "f00 = ", f00  # needs to be at q=0

#        E0=average(Ex)  # test
    E0=Ex

    fa=fa_asf(element,E0)
    sigh=debye_sears(element,T,qx) #;print "sqrt(dwf) = ", sigh   # square root of DWF i.e.  exp(-W)

    expF=expFh_dia(h,k,l)                    #; print expF
    Fh = sigh*(f0h+fa)*expF                  #; print Fh
    Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
    F0 = (f00+fa)*expFh_dia(0,0,0)           #; print F0   DWF=1 at q=0 ???

    V=a**3.0
    Chi0  = -re*F0/(pi*V)*(hpl*cl/E0)**2.0   #; print "chi_{0} = ",  Chi0
    Chih  = -re*Fh/(pi*V)*(hpl*cl/E0)**2.0   #; print "chi_{h} = ",  Chih
    Chih_ = -re*Fh_/(pi*V)*(hpl*cl/E0)**2.0  #; print "chi_{-h} = ", Chih_

##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)
    wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh
    # energy width for thick non-absorbing crystal (2.119)
    eps_s=4.0*re*dh**2.0/(pi*V)*abs(P*Fh)
    eps=eps_s/sqrt(abs(bh))

    K0_v=2.0*pi/(hpl*cl/Ex)
    alphah=H0/K0_v*(H0/K0_v-2.0*sin(thb))
    alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))

    #now !!!important!!! need to choose the sign of coren so that imag(eps1-eps2)>0
    coren=sqrt(alpha_pr**2.0+P**2.0*bh*Chih*Chih_) #; print "imag(coren) ", imag(coren)

    if imag(coren)<=0:
        coren=-coren

    eps1=Chi0-alpha_pr+coren #; print eps1
    eps2=Chi0-alpha_pr-coren #; print eps2
    kap1=0.5*eps1*K0/G0    #; print kap1
    kap2=0.5*eps2*K0/G0    #; print kap2
    R1=(eps1-Chi0)/(P*Chih_) #; print R1
    R2=(eps2-Chi0)/(P*Chih_) #; print R2

#        print "dc = ", dc
    Delta=(kap1-kap2)*dc        #; print dkap

    if Gh<0:
#          print "Bragg case"
        den0=R2-R1*exp(1.0j*Delta)  #; print den0
        t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
        r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h

#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    elif Gh>0:
#          print "Laue case"
        den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
        t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
        r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")

    else:
        print "Gh = 0, check input parameters"

    Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
    Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)

#        deny=2.0*abs(P*Chih)*sqrt(abs(bh))
#        y=alphah*bh+Chi0*(1.0-bh)/deny

    # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
    omeg0=2.0*sqrt(eps_s)
    bth=abs(0.5*pi-thb) #; print "bth", bth

    if average(bth) < average(omeg0):   #backscattering
        dth_s=omeg0
        dth=dth_s
    else:                               #non-backscattering
        dth_s=eps_s*tan(thb)
        dth=dth_s/sqrt(abs(bh))

    dth_pr=sqrt(abs(bh))*dth_s
    eps_pr=eps_s*sqrt(abs(bh))

    #de=sqrt(G0*abs(Gh))/(K0*abs(P*Chih)) # extinction length (2.90)

    #if rank(thb)==1 or rank(Ex)==1:         # check if vectors
       #de=de*imag(1.0/sqrt(y**2.0-1.0))     # extinction length (2.89)
       #de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    return [Chi0,Chih,Chih_,wh_s,wh,eps_s,eps,dth_s,dth,de,Tplot,Rplot,dth_pr,eps_pr]


######################################################################################################
# NEW FORM
######################################################################################################

def dtxrd0(thb,eta,phi,dc,Ex,P,crystalx):

    [a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element]=crystalx
###############################################################################
##  GEOMETRY
###############################################################################
    lamx=hpl*cl/Ex
    K0=2.0*pi/lamx   # ; print "K0 = ", K0
    H0=2.0*pi/dh     # ; print "H0 = ", H0
    G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
    Gh=G0-H0/K0*cos(eta)                            #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                    #; print "Gh = ", Gh # same actually
    bh=G0/Gh
###############################################################################
## STRUCTURE FACTOR and SUSCEPTIBILITIES
###############################################################################
########
    qx=0.5/dh #; print "qx (dh) = ", qx
#        qx1=sin(thb)/(hpl*cl/Ex); print "qx (lam) = ", qx1
    E0=Ex

    fa=fa_asf(element,E0)
    Fh = sigh*(f0h+fa)*expF                  #; print Fh
    Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
    F0 = (f00+fa)*expF0           #; print F0   DWF=1 at q=0 ???

    Chi0  = -re*F0/(pi*V)*(hpl*cl/E0)**2.0   #; print "chi_{0} = ",  Chi0
    Chih  = -re*Fh/(pi*V)*(hpl*cl/E0)**2.0   #; print "chi_{h} = ",  Chih
    Chih_ = -re*Fh_/(pi*V)*(hpl*cl/E0)**2.0  #; print "chi_{-h} = ", Chih_

##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)
    wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh
    # energy width for thick non-absorbing crystal (2.119)
    eps_s=4.0*re*dh**2.0/(pi*V)*abs(P*Fh)
    eps=eps_s/sqrt(abs(bh))

    K0_v=2.0*pi/(hpl*cl/Ex)
    alphah=H0/K0_v*(H0/K0_v-2.0*sin(thb))
    alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))

    #now !!!important!!! need to choose the sign of coren so that imag(eps1-eps2)>0
    coren=sqrt(alpha_pr**2.0+P**2.0*bh*Chih*Chih_) #; print "imag(coren) ", imag(coren)

    if imag(coren)<=0:
        coren=-coren

    eps1=Chi0-alpha_pr+coren #; print eps1
    eps2=Chi0-alpha_pr-coren #; print eps2
    kap1=0.5*eps1*K0/G0    #; print kap1
    kap2=0.5*eps2*K0/G0    #; print kap2
    R1=(eps1-Chi0)/(P*Chih_) #; print R1
    R2=(eps2-Chi0)/(P*Chih_) #; print R2

#        print "dc = ", dc
    Delta=(kap1-kap2)*dc        #; print dkap

    if Gh<0:
#          print "Bragg case"
        den0=R2-R1*exp(1.0j*Delta)  #; print den0
        t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
        r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h

#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    elif Gh>0:
#          print "Laue case"
        den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
        t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
        r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")

    else:
        print "Gh = 0, check input parameters"

    Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
    Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)

#        deny=2.0*abs(P*Chih)*sqrt(abs(bh))
#        y=alphah*bh+Chi0*(1.0-bh)/deny

    # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
    omeg0=2.0*sqrt(eps_s)
    bth=abs(0.5*pi-thb) #; print "bth", bth

    if average(bth) < average(omeg0):   #backscattering
        dth_s=omeg0
        dth=dth_s
    else:                               #non-backscattering
        dth_s=eps_s*tan(thb)
        dth=dth_s/sqrt(abs(bh))

    dth_pr=sqrt(abs(bh))*dth_s
    eps_pr=eps_s*sqrt(abs(bh))

    #de=sqrt(G0*abs(Gh))/(K0*abs(P*Chih)) # extinction length (2.90)

    #if rank(thb)==1 or rank(Ex)==1:         # check if vectors
       #de=de*imag(1.0/sqrt(y**2.0-1.0))     # extinction length (2.89)
       #de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    Chi=[Chi0,Chih,Chih_]
    return [Chi,wh_s,wh,eps_s,eps,dth_s,dth,de,Tplot,Rplot,dth_pr,eps_pr]


######################################################################################################
# NEW FORM with CHI function
# 2beam as of v0.26 and thruput4 (v0.14)
######################################################################################################

def dtxrd1(thb,eta,phi,dc,Ex,P,crystalx):
    [[Chi0,Chih,Chih_],dh]=crystalx
###############################################################################
##  GEOMETRY
###############################################################################
    Eb=0.5*hpl*cl/dh
    lamx=hpl*cl/Ex
    K0=2.0*pi/lamx   # ; print "K0 = ", K0
    H0=2.0*pi/dh     # ; print "H0 = ", H0
    G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
    Gh=G0-H0/K0*cos(eta)                            #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                    #; print "Gh = ", Gh # same actually
    bh=G0/Gh
##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)
    wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh
    # energy width for thick non-absorbing crystal (2.119)
    eps_s=4.0*dh**2.0/(lamx)**2.0*abs(P*Chih)
    eps=eps_s/sqrt(abs(bh))

    K0_v=2.0*pi/(hpl*cl/Ex)
    alphah=H0/K0_v*(H0/K0_v-2.0*sin(thb))
    alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))

    #now !!!important!!! need to choose the sign of coren so that imag(eps1-eps2)>0
    coren=sqrt(alpha_pr**2.0+P**2.0*bh*Chih*Chih_) #; print "imag(coren) ", imag(coren)

    if imag(coren)<=0:
        coren=-coren

    eps1=Chi0-alpha_pr+coren #; print eps1
    eps2=Chi0-alpha_pr-coren #; print eps2
    kap1=0.5*eps1*K0/G0    #; print kap1
    kap2=0.5*eps2*K0/G0    #; print kap2
    R1=(eps1-Chi0)/(P*Chih_) #; print R1
    R2=(eps2-Chi0)/(P*Chih_) #; print R2

#        print "dc = ", dc
    Delta=(kap1-kap2)*dc        #; print dkap

    if Gh<0:
#          print "Bragg case"
        den0=R2-R1*exp(1.0j*Delta)  #; print den0
        t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
        r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h

#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    elif Gh>0:
#          print "Laue case"
        den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
        t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
        r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")

    else:
        print "Gh = 0, check input parameters"

    Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
    Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)

#        deny=2.0*abs(P*Chih)*sqrt(abs(bh))
#        y=alphah*bh+Chi0*(1.0-bh)/deny

    # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
    omeg0=2.0*sqrt(eps_s)
    bth=abs(0.5*pi-thb) #; print "bth", bth

    if average(bth) < average(omeg0):   #backscattering
        eps_b=(Ex-Eb)/Eb               # (2.128-2.129)
        bthc2=2.0*(eps_b-wh_s)
        if bthc2 > eps_s:
            dth_s=sqrt(bthc2+eps_s)-sqrt(bthc2-eps_s)
        else:
            dth_s=2.0*sqrt(bthc2+eps_s)
        #
#             dth_s=omeg0               #simplified (2.130)
        dth=dth_s
    else:                               #non-backscattering
        dth_s=eps_s*tan(thb)
        dth=dth_s/sqrt(abs(bh))

    dth_pr=sqrt(abs(bh))*dth_s
    eps_pr=eps_s*sqrt(abs(bh))

    #de=sqrt(G0*abs(Gh))/(K0*abs(P*Chih)) # extinction length (2.90)

    #if rank(thb)==1 or rank(Ex)==1:         # check if vectors
       #de=de*imag(1.0/sqrt(y**2.0-1.0))     # extinction length (2.89)
       #de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    ep=[eps_s,eps,eps_pr]
    dt=[dth_s,dth,dth_pr]

#        return [wh_s,wh,eps_s,eps,dth_s,dth,de,Tplot,Rplot,dth_pr,eps_pr]
    return [wh_s,wh,ep,dt,de,Tplot,Rplot]
