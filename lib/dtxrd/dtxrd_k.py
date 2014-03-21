#!/usr/bin/env python

'''
a subroutine to calculate a set of parameters for transmitted/reflected a monochromatic wave

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

################################################################
# v 0.24 #######################################################
######################################################################################################################
# history:
# v0.24 02/21/2013  - implement 1beam case for energies above Eb
######################################################################################################################
from numpy import *
from fh import *
from constants import *

def dtxrd_k(k,eta,dc,Ex,P,crystalx):
    [kx,ky,kz]=k
    [a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element]=crystalx

###############################################################################
##  GEOMETRY
###############################################################################
    lamx=hpl*cl/Ex
    K0=2.0*pi/lamx    #; print "K0 = ", K0
    H0=2.0*pi/dh

    G0=kx*sin(eta)-kz*cos(eta)                  # ; print "G0 = ", G0
#        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
###############################################################################
## STRUCTURE FACTOR and SUSCEPTIBILITIES
###############################################################################
########
    fa=fa_asf(element,Ex)
    Fh = sigh*(f0h+fa)*expF                  #; print Fh
    Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
    F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???

    Chi0  = -re*F0/(pi*V)*lamx**2.0   #; print "chi_{0} = ",  Chi0
    Chih  = -re*Fh/(pi*V)*lamx**2.0   #; print "chi_{h} = ",  Chih
    Chih_ = -re*Fh_/(pi*V)*lamx**2.0  #; print "chi_{-h} = ", Chih_

    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)

    if Ex < Eb: #*(1.0+wh_s):   # 1 beam case  (p.58)
#        print " 1-beam case"

#       print "K0 = ", K0
#       print "Chi0 = ", Chi0
#       print "dc = ", dc
        kap0=0.5*K0*Chi0/abs(G0)  #abs(G0) needs to be used - look (2.52) ;print "kap0 = ", kap0
        kx_pr=kx; ky_pr=ky; kz_pr=kz #+kap0

        t00=exp(1.0j*kap0*dc)
        r0h=0.0
        Tplot=(abs(t00))**2.0    #; print "Tplot = ", Tplot
        Rplot=(abs(r0h))**2.0
        thc=arcsin(kz)
        dth=-1
        eps=-1
        bh=-1

    else:                 # 2 beam case

##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
####### direction cosine
        Gh=G0-H0/K0*cos(eta)                             #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                     #; print "Gh = ", Gh # same actually
        bh=G0/Gh #; print 'bh =', bh

        wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh
        thc=arcsin((Eb/Ex)*(1.0+wh))       # optimal angle

        # energy width for thick non-absorbing crystal (2.119)
        eps_s=4.0*re*dh**2.0/(pi*V)*abs(P*Fh)
        eps=eps_s/sqrt(abs(bh))

#        K0_v=2.0*pi/(hpl*cl/Ex)
        alphah=H0/K0*(H0/K0+2.0*kz)

####### Momentum transfer (2.17) ######################
        coren0=sign(Gh)*sqrt(Gh**2.0-alphah)
        delh=K0*(-Gh+coren0) # ; print "delh = ", delh          # Momentum transfer
##########################################################################################
######### Analytical solution for 2-beam case ############################################
        alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))  # (2.60)
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

#        coren0=sqrt(Gh**2.0-alphah)

        if Gh<0:
#          print "Bragg case"
            den0=R2-R1*exp(1.0j*Delta)  #; print den0
            t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
            r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h

#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        elif Gh>0:
#          print "Laue case"
            den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
            t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
            r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        else:
            print "Gh = 0, check input parameters"

        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)

######### Exit wavevector  ##########################################################
        kx_pr=kx+sin(eta)*delh/K0
        ky_pr=ky
        kz_pr=kz+H0/K0-cos(eta)*delh/K0
######### input/exit parameters #####################################################
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thc) #; print "bth", bth

        if average(bth) < average(omeg0):   #backscattering
            dth_s=omeg0
            dth=dth_s
        else:                               #non-backscattering
            dth_s=eps_s*tan(thc)
            dth=dth_s/sqrt(abs(bh))

        dth_pr=sqrt(abs(bh))*dth_s
        eps_pr=eps_s*sqrt(abs(bh))

        #de=sqrt(G0*abs(Gh))/(K0*abs(P*Chih)) # extinction length (2.90)
        de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    return [[kx_pr,ky_pr,kz_pr],thc,dth,eps,bh,Tplot,Rplot]



################################################################################
## Calculates Chi and everything
################################################################################

def dtxrd0_k(k,eta,dc,Ex,P,crystalx):
    [kx,ky,kz]=k
    [a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element]=crystalx

###############################################################################
##  GEOMETRY
###############################################################################
    lamx=hpl*cl/Ex
    K0=2.0*pi/lamx    #; print "K0 = ", K0
    H0=2.0*pi/dh

    G0=kx*sin(eta)-kz*cos(eta)                  # ; print "G0 = ", G0
#        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
###############################################################################
## STRUCTURE FACTOR and SUSCEPTIBILITIES
###############################################################################
########
    fa=fa_asf(element,Ex)
    Fh = sigh*(f0h+fa)*expF                  #; print Fh
    Fh_= sigh*(f0h+fa)*expF.conjugate()      #; print Fh_
    F0 = (f00+fa)*expF0                      #; print F0   DWF=1 at q=0 ???

    Chi0  = -re*F0/(pi*V)*lamx**2.0   #; print "chi_{0} = ",  Chi0
    Chih  = -re*Fh/(pi*V)*lamx**2.0   #; print "chi_{h} = ",  Chih
    Chih_ = -re*Fh_/(pi*V)*lamx**2.0  #; print "chi_{-h} = ", Chih_

    Chi=[Chi0,Chih,Chih_]

    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)

    if Ex < Eb: #*(1.0+wh_s):   # 1 beam case  (p.58)
#        print " 1-beam case"

#       print "K0 = ", K0
#       print "Chi0 = ", Chi0
#       print "dc = ", dc
        kap0=0.5*K0*Chi0/abs(G0)  #abs(G0) needs to be used - look (2.52) ;print "kap0 = ", kap0
        kx_pr=kx; ky_pr=ky; kz_pr=kz #+kap0

        t00=exp(1.0j*kap0*dc)
        r0h=0.0
        Tplot=(abs(t00))**2.0    #; print "Tplot = ", Tplot
        Rplot=(abs(r0h))**2.0
        thc=arcsin(kz)
        dth=-1
        eps=-1
        bh=-1

    else:                 # 2 beam case

##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
####### direction cosine
        Gh=G0-H0/K0*cos(eta)                             #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                     #; print "Gh = ", Gh # same actually
        bh=G0/Gh #; print 'bh =', bh

        wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh
        thc=arcsin((Eb/Ex)*(1.0+wh))       # optimal angle

        # energy width for thick non-absorbing crystal (2.119)
        eps_s=4.0*re*dh**2.0/(pi*V)*abs(P*Fh)
        eps=eps_s/sqrt(abs(bh))

#        K0_v=2.0*pi/(hpl*cl/Ex)
        alphah=H0/K0*(H0/K0+2.0*kz)

####### Momentum transfer (2.17) ######################
        coren0=sign(Gh)*sqrt(Gh**2.0-alphah)
        delh=K0*(-Gh+coren0) # ; print "delh = ", delh          # Momentum transfer
##########################################################################################
######### Analytical solution for 2-beam case ############################################
        alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))  # (2.60)
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

#        coren0=sqrt(Gh**2.0-alphah)

        if Gh<0:
#          print "Bragg case"
            den0=R2-R1*exp(1.0j*Delta)  #; print den0
            t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
            r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h

#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        elif Gh>0:
#          print "Laue case"
            den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
            t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
            r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        else:
            print "Gh = 0, check input parameters"

        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)

######### Exit wavevector  ##########################################################
        kx_pr=kx+sin(eta)*delh/K0
        ky_pr=ky
        kz_pr=kz+H0/K0-cos(eta)*delh/K0
######### input/exit parameters #####################################################
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thc) #; print "bth", bth

        if average(bth) < average(omeg0):   #backscattering
            dth_s=omeg0
            dth=dth_s
        else:                               #non-backscattering
            dth_s=eps_s*tan(thc)
            dth=dth_s/sqrt(abs(bh))

        dth_pr=sqrt(abs(bh))*dth_s
        eps_pr=eps_s*sqrt(abs(bh))

        #de=sqrt(G0*abs(Gh))/(K0*abs(P*Chih)) # extinction length (2.90)
        de=1.0/(imag(kap1-kap2))              # extinction length ("precise")

    return [[kx_pr,ky_pr,kz_pr],thc,dth,eps,bh,Tplot,Rplot,Chi]


###############################################################################
## Fast - calculates only k_out, T and R
###############################################################################

def dtxrd1_k(k,eta,dc,Ex,P,crystalx,Chi):
    [kx,ky,kz]=k
    [a,dh,Eb,f0h,f00,expF0,expF,sigh,V,element]=crystalx
    [Chi0,Chih,Chih_]=Chi
###############################################################################
##  GEOMETRY
###############################################################################
#    lamx=hpl*cl/Ex
#    K0=2.0*pi/lamx    #; print "K0 = ", K0
    K0=2.0*pi*Ex/(hpl*cl)
    H0=2.0*pi/dh
    G0=kx*sin(eta)-kz*cos(eta)                  # ; print "G0 = ", G0
#
#    wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)

    if Ex < Eb: #*(1.0+wh_s):   # 1 beam case  (p.58)
        print " 1-beam case"
        kap0=0.5*K0*Chi0/abs(G0)  #abs(G0) needs to be used - look (2.52) ;print "kap0 = ", kap0
        kx_pr=kx; ky_pr=ky; kz_pr=kz #+kap0
        t00=exp(1.0j*kap0*dc)
        r0h=0.0
        Tplot=(abs(t00))**2.0    #; print "Tplot = ", Tplot
        Rplot=(abs(r0h))**2.0
        thc=arcsin(kz)

    else:                 # 2 beam case

##########################################################################################
####### FOR DYNAMICAL BRAGGs LAW
##########################################################################################
####### direction cosine
        Gh=G0-H0/K0*cos(eta)                             #; print "Gh = ", Gh
        bh=G0/Gh #; print 'bh =', bh
        alphah=H0/K0*(H0/K0+2.0*kz)
####### Momentum transfer (2.17) ######################
        coren0=sign(Gh)*sqrt(Gh**2.0-alphah)
        delh=K0*(-Gh+coren0) # ; print "delh = ", delh          # Momentum transfer
##########################################################################################
######### Analytical solution for 2-beam case ############################################
        alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))  # (2.60)
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
        Delta=(kap1-kap2)*dc        #; print dkap
#        coren0=sqrt(Gh**2.0-alphah)

        if Gh<0: # print "Bragg case"
            den0=R2-R1*exp(1.0j*Delta)  #; print den0
            t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00
            r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0  #; print r0h
#        den0=R2*exp(-1.0j*Delta)-R1  #; print den0
#        t00=exp(1.0j*kap2*dc)*(R2-R1)/den0             #; print t00
#        r0h=R1*R2*(exp(-1.0j*Delta)-1.0)/den0  #; print r0h
        elif Gh>0: # print "Laue case"
            den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
            t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
            r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
        else:
            print "Gh = 0, check input parameters"
#
        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)
######### Exit wavevector  ###################################################################
        kx_pr=kx+sin(eta)*delh/K0
        ky_pr=ky
        kz_pr=kz+H0/K0-cos(eta)*delh/K0

    return [[kx_pr,ky_pr,kz_pr],Tplot,Rplot]
