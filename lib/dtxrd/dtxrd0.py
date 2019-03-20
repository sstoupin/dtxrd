#!/usr/bin/env python

'''
a subroutine to calculate a set of parameters for transmitted/reflected a monochromatic wave

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

from numpy import *
from fh import *
from constants import *
from deriv import *
from numpy.lib import scimath as SM
#-----------------------------------------------------------------------------------------------
# v 0.20
#-----------------------------------------------------------------------------------------------
# v 0.2 added a subroutine dtxtd1_D0Dh to calculate radiation field inside the crystal 
#       rename this main subroutine to dtxrd0.py to avoid stuff like dtxrd.dtxrd
#-----------------------------------------------------------------------------------------------
#

def dtxrd(element,h,k,l,thb,eta,phi,a,dh,T,dc,Ex,P):

###############################################################################
##  GEOMETRY
############################################################################### 
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
          print("Gh = 0, check input parameters")
        
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
        else:             		    #non-backscattering
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
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
          print("Gh = 0, check input parameters")   
        
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
        else:             		    #non-backscattering
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
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx 
###############################################################################
##  GEOMETRY
############################################################################### 
        Eb=0.5*hpl*cl/dh
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih) #*Chih_)  # Authier p. 102 - same result ! 01/25/2015

        elif Gh>0:
#          print "Laue case"
          den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
          t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
          r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
          de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih)   # Authier p. 102 - same result!      01/25/2015
          
        else:
          print("Gh = 0, check input parameters")   
        
        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)
        
#        deny=2.0*abs(P*Chih)*sqrt(abs(bh))        
#        y=alphah*bh+Chi0*(1.0-bh)/deny        
                
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thb) #; print "bth", bth
        
        if average(bth) < average(omeg0):   #backscattering
             # calculate Er: ##########################################################################################################################
#             eta=eta+1.0e-6
#             denr=1.0+sqrt(1.0+4.0*wh_s*(tan(eta))**2.0)
#             bthr=abs(arctan(2.0*wh_s*tan(eta)/denr))
#             lamr=2.0*dh*cos(bthr)*(1.0-tan(bthr)/tan(eta))
#             Er=hpl*cl/lamr
             eps_b=(Ex-Eb)/Eb         # (2.128-2.129)                                      
#             eps_b=(Ex-Er)/Er         # better approximation  (Ex-Er)/Er - attempted - actually gives different values!!!
                                       # the RuntimeWarning is due to calculation of eps_b at different energies beyond the Darwin table
                                       # this message will be suppressed by taking the abs. value |bthc2+eps_s| - does not affect data within 
                                       # the Darwin table
             ########################################################################################################################################## 
             bthc2=2.0*(eps_b-wh_s)
             if bthc2 > eps_s:
                dth_s=sqrt(bthc2+eps_s)-sqrt(bthc2-eps_s)
             else:             
                dth_s=2.0*sqrt(abs(bthc2+eps_s))
             ###########################################################################################################################################
#             dth_s=omeg0               #simplified (2.130)
             dth=dth_s
        else:             		    #non-backscattering
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
########################################################################################################################################################
########################################################################################################################################################
# zero absorption SOLUTION
########################################################################################################################################################
def dtxrd1z(thb,eta,phi,dc,Ex,P,crystalx):
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx 
        Chi0=real(Chi0)
        Chih=real(Chih)
        Chih_=real(Chih_)
###############################################################################
##  GEOMETRY
############################################################################### 
        Eb=0.5*hpl*cl/dh
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
        coren=SM.sqrt(alpha_pr**2.0+P**2.0*bh*Chih*Chih_) #; print "coren", coren

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
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih) #*Chih_)  # Authier p. 102 - same result ! 01/25/2015

        elif Gh>0:
#          print "Laue case"
          den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
          t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
          r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
          de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih)   # Authier p. 102 - same result!      01/25/2015
          
        else:
          print("Gh = 0, check input parameters")
        
        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)
        
#        deny=2.0*abs(P*Chih)*sqrt(abs(bh))        
#        y=alphah*bh+Chi0*(1.0-bh)/deny        
                
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thb) #; print "bth", bth
        
        if average(bth) < average(omeg0):   #backscattering
             # calculate Er: ##########################################################################################################################
#             eta=eta+1.0e-6
#             denr=1.0+sqrt(1.0+4.0*wh_s*(tan(eta))**2.0)
#             bthr=abs(arctan(2.0*wh_s*tan(eta)/denr))
#             lamr=2.0*dh*cos(bthr)*(1.0-tan(bthr)/tan(eta))
#             Er=hpl*cl/lamr
             eps_b=(Ex-Eb)/Eb         # (2.128-2.129)                                      
#             eps_b=(Ex-Er)/Er         # better approximation  (Ex-Er)/Er - attempted - actually gives different values!!!
                                       # the RuntimeWarning is due to calculation of eps_b at different energies beyond the Darwin table
                                       # this message will be suppressed by taking the abs. value |bthc2+eps_s| - does not affect data within 
                                       # the Darwin table
             ########################################################################################################################################## 
             bthc2=2.0*(eps_b-wh_s)
             if bthc2 > eps_s:
                dth_s=sqrt(bthc2+eps_s)-sqrt(bthc2-eps_s)
             else:             
                dth_s=2.0*sqrt(abs(bthc2+eps_s))
             ###########################################################################################################################################
#             dth_s=omeg0               #simplified (2.130)
             dth=dth_s
        else:             		    #non-backscattering
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
##############################################################################################################################################################
# Reflectivity/Transmissivity solutions of Darwin-Hamilton equations for mosaic crystals
# needs 1d array for sig_v
##############################################################################################################################################################
def darwin_hamilton(thc,eta,phi,dc,Ec,P,crystalx,sig_v):
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx
        G0=cos(thc)*sin(eta)*cos(phi)+sin(thc)*cos(eta)
        lamx=hpl*cl/Ec
        K0 = 2.0*pi/lamx
        H0 = 2.0*pi/dh
        Gh=G0-H0/K0*cos(eta)
        #                 
        delta = -0.5*real(Chi0)
        beta=0.5*imag(Chi0)
        mu0 = 2.0*beta*K0  #[Angstrom] 
        #-------------------------------------------------------
        # criteria parameters for primary/secondary extinction:
        le = lamx*sqrt(G0*abs(Gh))/(abs(P)*sqrt(Chih*Chih_))
        le = real(le)
        #A = pi*K0*P*abs(Chih)*t0/sqrt(abs(G0*Gh)) # this is crystal block size divided by the extinction length in perfect crystal
        A  = dc/le
        Gs_mu = max(sig_v)/mu0   # simply scattering to absorption ratio (in the same direction of propagation)
        #-------------------------------------------------------
        Gp_=0.5*(1.0/G0 + 1.0/abs(Gh))
        Gm_=0.5*(1.0/G0 - 1.0/abs(Gh))
        #eps2 = sqrt(Gm_**2.0*(mu0+sig_v)**2.0 + sig_v**2.0/(G0*abs(Gh)))
        #epsb = sqrt(Gp_**2.0*(mu0+sig_v)**2.0 - sig_v**2.0/(G0*abs(Gh)))
        
        if Gh < 0:  # Bragg case
            epsb = sqrt(Gp_**2.0*(mu0+sig_v)**2.0 - sig_v**2.0/(G0*abs(Gh)))
            Rmos = sig_v/G0*sinh(epsb*dc)/(Gp_*(mu0+sig_v)*sinh(epsb*dc) + epsb*cosh(epsb*dc))
            Tmos = epsb*exp(-Gm_*(mu0+sig_v)*dc)/(Gp_*(mu0+sig_v)*sinh(epsb*dc) + epsb*cosh(epsb*dc))
            opt_d = 'NA'
            Tpmos = exp(-(mu0+sig_v)*dc/G0) # transm. for perfect mosaic crystal
            Rpmos = sig_v/G0/(mu0*Gp_)*sinh(mu0*Gp_*dc)/exp(mu0*Gp_*dc) # reflectivity for perfect mosaic crystal
            #
        elif Gh > 0: #Laue case
            eps2 = sqrt(Gm_**2.0*(mu0+sig_v)**2.0 + sig_v**2.0/(G0*abs(Gh)))
            #
            Rmos = sig_v/(G0*eps2)*exp(-Gp_*(mu0+sig_v)*dc)*sinh(eps2*dc)
            Tmos = 1.0/eps2*exp(-Gp_*(mu0+sig_v)*dc)*(-(mu0+sig_v)*Gm_*sinh(eps2*dc) + eps2*cosh(eps2*dc))
            eps2_max = sqrt(Gm_**2.0*(mu0+max(sig_v))**2.0 + max(sig_v)**2.0/(G0*abs(Gh)))
            opt_d = 1.0/eps2_max*arctanh(eps2_max/(Gp_*(mu0+max(sig_v))))
            Tpmos = exp(-(mu0+sig_v)*dc/G0)                           # maybe not very precise in the first order on sig_v
            Rpmos = exp(-mu0*dc*Gp_)*sig_v/(Gm_*mu0)*sinh(Gm_*mu0*dc) # reflectivity for perfect mosaic crystal (asym.)
            #
        return [A,Gs_mu,opt_d,Tpmos,Rpmos,Tmos,Rmos]        
##############################################################################################################################################################
##############################################################################################################################################################
def bragg_laue(thb,eta,phi,dc,Ex,P,crystalx,t):  # t- distance from entrance point to lateral crystal face
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx 
###############################################################################
##  GEOMETRY
############################################################################### 
        Eb=0.5*hpl*cl/dh
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
#        Gh=G0-2.0*sin(thb)*cos(eta)                    #; print "Gh = ", Gh # same actually
        bh=G0/Gh
        #
        delta = -0.5*real(Chi0)
        beta=0.5*imag(Chi0)
        mu0 = 2.0*beta*K0  #[Angstrom]         
        
        # calculate deviation parameter
        K0_v=2.0*pi/(hpl*cl/Ex)
        alphah=H0/K0_v*(H0/K0_v-2.0*sin(thb))                                
        alpha_pr=0.5*(alphah*bh+Chi0*(1.0-bh))
        
        #---------------------------------------------------------------
                                
        #wh_s=-2.0*real(Chi0)*(dh/lamx)**2.0  #; print "wh(s) = ", wh_s #wh_s=average(wh_s)        
        #wh=0.5*wh_s*(bh-1.0)/bh              #; print "wh    = ", wh                
        # energy width for thick non-absorbing crystal (2.119)
        #eps_s=4.0*dh**2.0/(lamx)**2.0*abs(P*Chih)
        #eps=eps_s/sqrt(abs(bh))                        
                
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
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih) #*Chih_)  # Authier p. 102 - same result ! 01/25/2015

        elif Gh>0:
#          print "Laue case"
          den0=R1-R2
#          t00=exp(1.0j*kap2*dc)*(R1-R2*exp(1.0j*Delta))/den0
#          r0h=R1*R2*exp(1.0j*kap2*dc)*(1.0-exp(1.0j*Delta))/den0
          t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
          r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
          de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")
#          de=0.5/pi*lamx*sqrt(G0*abs(Gh))/abs(Chih)   # Authier p. 102 - same result!      01/25/2015
          
        else:
          print("Gh = 0, check input parameters")   
        
        
        # reduced deviation parameter:        
        deny=2.0*abs(P*Chih)*sqrt(abs(bh))        
        y=alphah*bh+Chi0*(1.0-bh)/deny
        # Authier 5.46 #################################################
        lam0 = lamx*sqrt(G0*abs(Gh))/(P*sqrt(Chih*Chih_))
        mu_e = mu0 + 2.0*pi*G0*imag((y + sqrt(y**2.0 - 1.0))/lam0)
        
        T_a=(abs(t00))**2.0          #t00*t00.conjugate()
        R_a = (abs(r0h))**2.0/abs(bh)  # refl. on entrance face
        angle = (thb+eta)/r2d
        R_d = arctan(angle)*Gh*abs(r0h)**2.0*exp(-mu_e*t/G0) #????
        T_d = arctan(angle)*G0/sin(thb/r2d+eta/r2d)*exp(-mu_e*t/G0)  #????
        #                                
        return [T_a,R_a,T_d,R_d] 
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
# RETURNS FIELD inside the crystal along z
# from z=0 to z=d in ns+1 steps
###############################################################################
def dtxrd1_D0Dh(thb,eta,phi,dc,Ex,P,crystalx,zv):
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx 
###############################################################################
##  GEOMETRY
############################################################################### 
        Eb=0.5*hpl*cl/dh
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
        
        if Gh<0:
#          print "Bragg case"
          Delta=(kap1-kap2)*dc                         #; print dkap
          den0=R2-R1*exp(1.0j*Delta)                   #; print den0
          t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00        
          r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0         #; print r0h                   
          de=1.0/(imag(kap1-kap2))                     # extinction length ("precise")
          #-------------------------------------------------------------------------------- 
          # Radiation fields inside the crystal:
          #--------------------------------------------------------------------------------
          denD = R1*exp(-1.0j*kap2*dc)-R2*exp(-1.0j*kap1*dc)
          D0 = (R1*exp(1.0j*kap2*(zv-dc))-R2*exp(1.0j*kap1*(zv-dc)))/denD
          Dh = R1*R2*(exp(1.0j*kap2*(zv-dc))-exp(1.0j*kap1*(zv-dc)))/denD
          #
        elif Gh>0:
#          print "Laue case"
          den0=R1-R2
          t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
          r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
          de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")
          #-------------------------------------------------------------------------------- 
          # Radiation fields inside the crystal:
          #--------------------------------------------------------------------------------
          denD = R1-R2
          D0 = (R1*exp(1.0j*kap2*zv)-R2*exp(1.0j*kap1*zv))/denD
          Dh = R1*R2*(exp(1.0j*kap2*zv)-exp(1.0j*kap1*zv))/denD
          #
        else:
          print("Gh = 0, check input parameters")
        
        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)
        #
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thb) #; print "bth", bth
        
        if average(bth) < average(omeg0):   #backscattering
             # calculate Er: ##########################################################################################################################
#             eta=eta+1.0e-6
#             denr=1.0+sqrt(1.0+4.0*wh_s*(tan(eta))**2.0)
#             bthr=abs(arctan(2.0*wh_s*tan(eta)/denr))
#             lamr=2.0*dh*cos(bthr)*(1.0-tan(bthr)/tan(eta))
#             Er=hpl*cl/lamr
             eps_b=(Ex-Eb)/Eb          # (2.128-2.129)                                      
#             eps_b=(Ex-Er)/Er          # better approximation 
             ##########################################################################################################################################
             bthc2=2.0*(eps_b-wh_s)
             if bthc2 > eps_s:
                dth_s=sqrt(bthc2+eps_s)-sqrt(bthc2-eps_s)
             else:
                dth_s=2.0*sqrt(abs(bthc2+eps_s)) # abs to suppress the RuntimeWarning at Ex outside the Darwin table             
             #
              #dth_s=omeg0                  #simplified (2.130)
             dth=dth_s
        else:             		    #non-backscattering
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
        GG=[G0,Gh]
        DD=[D0,Dh]
        #                                                      
        return [[wh_s,wh,ep,dt,de,Tplot,Rplot],[GG,DD]]

###############################################################################
# Calculates normalized electron yield by integrating along z
###############################################################################
def dtxrd1_sy(thb,eta,phi,dc,Ex,P,crystalx,Ls):
        [[Chi0,Chih,Chih_],dh,lp,ang]=crystalx 
        #######################################################################
        # Z vector for numerical calculation of the secondary yield
        ####################################################################### 
        if Ls < dc/3.0:
            zran=3.0*Ls 
        else: 
            zran=dc
        zstep=zran/float(100)
        zv = arange(0.0,zran,zstep) 
###############################################################################
##  GEOMETRY
############################################################################### 
        Eb=0.5*hpl*cl/dh
        lamx=hpl*cl/Ex
        K0=2.0*pi/lamx   # ; print "K0 = ", K0
        H0=2.0*pi/dh     # ; print "H0 = ", H0
        G0=cos(thb)*sin(eta)*cos(phi)+sin(thb)*cos(eta) #; print "G0 = ", G0
        Gh=G0-H0/K0*cos(eta)			        #; print "Gh = ", Gh
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
        
        if Gh<0:
#          print "Bragg case"
          Delta=(kap1-kap2)*dc                         #; print dkap
          den0=R2-R1*exp(1.0j*Delta)                   #; print den0
          t00=exp(1.0j*kap1*dc)*(R2-R1)/den0           #; print t00        
          r0h=R1*R2*(1.0-exp(1.0j*Delta))/den0         #; print r0h                   
          de=1.0/(imag(kap1-kap2))                     # extinction length ("precise")
          #-------------------------------------------------------------------------------- 
          # Radiation fields inside the crystal:
          #--------------------------------------------------------------------------------
          denD = R1*exp(-1.0j*kap2*dc)-R2*exp(-1.0j*kap1*dc)
          D0 = (R1*exp(1.0j*kap2*(zv-dc))-R2*exp(1.0j*kap1*(zv-dc)))/denD
          Dh = R1*R2*(exp(1.0j*kap2*(zv-dc))-exp(1.0j*kap1*(zv-dc)))/denD          
          #
          #----------------------------------------------------------------------------------
        elif Gh>0:
#          print "Laue case"
          den0=R1-R2
          t00=(R1*exp(1.0j*kap2*dc)-R2*exp(1.0j*kap1*dc))/den0
          r0h=R1*R2*(exp(1.0j*kap2*dc)-exp(1.0j*kap1*dc))/den0
          de=1.0/(abs(real(kap1-kap2)))              # extinction length ("precise")
          #-------------------------------------------------------------------------------- 
          # Radiation fields inside the crystal:
          #--------------------------------------------------------------------------------
          denD = R1-R2
          D0 = (R1*exp(1.0j*kap2*zv)-R2*exp(1.0j*kap1*zv))/denD
          Dh = R1*R2*(exp(1.0j*kap2*zv)-exp(1.0j*kap1*zv))/denD
          #
        else:
          print("Gh = 0, check input parameters")   
        ##################################################################################
        # REFLECTIVITY AND TRANSMISSIVITY
        ##################################################################################
        Tplot=(abs(t00))**2.0          #t00*t00.conjugate()
        Rplot=(abs(r0h))**2.0/abs(bh)  #;print 'R = ', Rplot #r0h*r0h.conjugate()/abs(bh)
        ##################################################################################
        # SECONDARY YIELD
        ##################################################################################        
        #
        D02=abs(D0)**2.0; [zz,dD02]=deriv(zv,D02)
        Dh2=abs(Dh)**2.0; [zz,dDh2]=deriv(zv,Dh2)
        kappaz=-G0*dD02-Gh*dDh2   # -Gh gets positive when Gh<0 (Bragg) and negative when Gh>0 (Laue)
        Pz=1.0*exp(-zz/Ls)        # for simplicity assume gain factor = 1.0
        syieldz=array(kappaz)*array(Pz)
        SYplot=zstep*sum(syieldz)                                                                
        #################################################################################
        # ANGULAR AND ENERGY ENTRANCE AND EXIT 
        #################################################################################
        # angular width for thick non-absorbing crystal (far from bc) (2.125-2.126)
        omeg0=2.0*sqrt(eps_s)
        bth=abs(0.5*pi-thb) #; print "bth", bth
        #        
        if average(bth) < average(omeg0):   #backscattering
             eps_b=(Ex-Eb)/Eb               # (2.128-2.129)             
             bthc2=2.0*(eps_b-wh_s)
             if bthc2 > eps_s:
                dth_s=sqrt(bthc2+eps_s)-sqrt(bthc2-eps_s)
             else:
                dth_s=2.0*sqrt(abs(bthc2+eps_s))  # abs - to suppress the RuntimeWarning at Ex outside of the Darwin table
             #
              #dth_s=omeg0                  #simplified (2.130)
             dth=dth_s
        else:             		    #non-backscattering
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
        GG=[G0,Gh]
        DD=[D0,Dh]
        #                                                      
        return [wh_s,wh,ep,dt,de,Tplot,Rplot,SYplot]
                                                        