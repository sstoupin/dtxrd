#!/usr/bin/env python

'''
x-ray diffraction calculator (2-beam case, perfect crystals)

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys                                                                      #@UnusedImport
#from numpy import *
#from scipy import *
from pylab import *                                                     #@UnusedWildImport

import os
if os.path.abspath(os.path.dirname(__file__)).split(os.sep)[-2] == 'lib':
    '''when running this script from the source directory'''
    sys.path.insert(0, os.path.abspath('..'))

from dtxrd.myio import coljoin, writeFile       #@UnusedImport
from dtxrd.curvestat import *                           #@UnusedWildImport
from dtxrd.thfind import *                                      #@UnusedWildImport
from dtxrd.dtxrd0 import *                                      #@UnusedWildImport
from dtxrd.constants import *                           #@UnusedWildImport
from dtxrd.chi import *                                         #@UnusedWildImport

__version__ = '0.29'
#-------------------------------------------------------------------------------------------
# SUMMARY OF CHANGES:
#-------------------------------------------------------------------------------------------
# v 0.11 help page on USAGE extended
# v 0.12 - change name to 2beam
# try to add 2 beam simulation
# v 0.13 basic "engine" for 2beam front surface incidence implemented
# v 0.14 DWF added
# v 0.15
# * anomalous str.f replaced interpolation with spline (see fh.py)
# * compare with Yuri's code - satisfactory so far
# *  spikes ??? divergence for Si? v0.15 - solved by taking f_pp with negative sign (see fh.py)
# v 0.16
# * centered ranges in accord. with dynamical Bragg law
# * added crystal thickness to input parameters v0.16 done
# v 0.161 added explicit path to data files (also in fh.py)
# v 0.17 added interpolation of Yuri's tables for f1 and f2 instead of cromer's values at fixed energies
# v 0.18 added energy width for thick non-absorbing crystal
# v 0.19 introduced subroutine dtxrd
# v 0.20 fixed problem with choosing the sign for "coren" in dtxrd
# v 0.21 dth and dE calculated using curvestat/convolution
# v 0.22 added Laue case to dtxrd
# v 0.23 write data to a file (Tplot, Rplot, rc0 (convolution))
# v 0.24 added polarization option
# v 0.25 now using myio stuff to write data
# v 0.26 new functions + include Al2O3
# v 0.27 included calculations of thr and Er
# v 0.28 included warning when expF is small. convolution with res. function is now an option -c
# v 0.29 better understood intrinsic angular dispersion rate - use (2.138) to calculate
#        more precise formula for the angular width in backscattering
#-------------------------------------------------------------------------------------------------------------
#  things to do
# 1. draw a figure of crystal and waves

def fatalError(msg):
    sys.stderr.write('Error: ')
    sys.stderr.write(str(msg))
    sys.stderr.write('\n')
    sys.exit(1)

def fatalIOError(err):
    if issubclass(err.__class__, IOError) and err.strerror and err.filename:
        err = '%s: %s' % (err.strerror, err.filename)
    fatalError(err)

def ParseArguments(args):
    try:
        from optik import OptionParser
    except ImportError:
        try:
            from optparse import OptionParser
        except ImportError:
            fatalError(
'This program requires Optik, available from http://optik.sourceforge.net/\n')

    USAGE = '%prog [OPTIONS...] element h k l eta phi T d flag theta(or Ex) \n'+'where:\n' \
    +'crystal - C (diamond), Si (silicon), Ge (germanium) or Al2O3 (sapphire)\n' \
    +'h k l - miller indices\n' \
    +'eta  - asymmetry angle \n' \
    +'phi  - azimuthal angle of incidence\n' \
    +'T -temperature [K] \n' \
    +'d - crystal thickness [mm] \n'\
    +'flag = a - calculate at given angle of incidence\n' \
    +'flag = e - calculate at given x-ray energy\n' \
    +'theta - glancing angle of incidence (flag = a) or Ex - x-ray energy [keV] (flag = e)'
    VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>\n' \
    +'the program calculates Bragg angle, input: material, reflection, temperature, energy'
    parser = OptionParser(usage=USAGE, version=VERSION)
#----------------------------------------------------------------------------------------------
    parser.add_option('-o', '--output',
            action='store',
            dest='output',
            default=None,
            help='write results to file F (defaults to stdout)',
            metavar='F')
    parser.add_option('-w', '--write',
            action='store',
            dest='write',
            default=None,
            help='write data to file (default do not write)',
            metavar='D')
    parser.add_option('-p', '--pi',
            action='store_const',
            const=1,
            dest='polarization',
            default=0,
            help='-p - pi polarization ')
    parser.add_option('-c', '--conv',
            action='store_const',
            const=1,
            dest='conv',
            default=0,
            help='-c - convolve with instrumental res. function')
#-----------------------------------------------------------------------------------------------
    opts, args = parser.parse_args(args)
    if len(args) != 10:
        parser.print_usage()
        sys.exit(1)

    try:
        args[9] = float(args[9])
    except ValueError:
        fatalError('Ex(keV) or theta_b (deg.)  must be a valid number')

    return opts, args
#-----------------------------------------------------------------------------------------------


def main():
    opts, args = ParseArguments(sys.argv[1:])

    if opts.output is not None:
        try:
            outFile = open(opts.output, 'w')
        except IOError, e:
            fatalIOError(e)
    else:
        outFile = sys.stdout

########
    h=float(args[1])
    k=float(args[2])
    l=float(args[3])
    eta=float(args[4])/r2d  #; print "cos(eta) =", cos(eta)
    phi=float(args[5])/r2d
    T=float(args[6])
    dc=1.0e7*float(args[7]) # [A] crystal thickness
    flag=args[8]
######## polarization factor
    P=1.0   # default value (sigma) - if pi - controlled by option -p
#########################################################################################################
## LATTICE PARAMETER AND BRAGG ENERGY
#########################################################################################################
    crystal,flagFh=chi(args[0],h,k,l,T,5.0e4)
    # break if structure amplitude is too small:
    if flagFh==0:
        fatalError('forbidden reflection: structure amplitude |exp(iHr)| for the chosen set of Miller indicies is too small (< 1e-6)')
    #------------------------------------------------------------------------------------
    dh=crystal[1]
    Eb=0.5*hpl*cl/dh
    thr,Er=thr_find(args[0],h,k,l,eta+1.0e-14,phi,T,dc,Eb)
######  BEGIN OUTPUT ###############################################################################
    outFile.write('##############################################################\n')
    outFile.write('##### 2BEAM v'+__version__+' ############################################\n')
    outFile.write('##### Author: Stanislav Stoupin ## sstoupin@aps.anl.gov ######\n')
    outFile.write('##############################################################\n')
    outFile.write('PARAMETERS OF THE H K L REFLECTION :\n')
#        outFile.write('a [A] = '+str(a)+' lattice parameter\n')
    outFile.write('d [A] = '+str(dh)+' d-spacing\n')
    outFile.write('Eb [keV]  = '+str(1.0e-3*Eb)+' Bragg Energy\n')
    outFile.write('thr [deg] = '+str(thr*r2d)+' Angle of the exact backscattering\n')
    outFile.write('Er [keV] = '+str(1.0e-3*Er)+' Energy of the exact backscattering\n')
    outFile.write('--------------------------------------------------------------\n')
####################################################################################################
    f1=plt.figure()
    if flag=='a':
        thb=float(args[9])/r2d
        thc=thb        # assume input angle as true central angle
        Ex=Eb/sin(thb) # approximate energy from kinematic Braggs Law

        if opts.polarization==1:
            P=cos(2.0*thc)

#             print "Ex = ", Ex
        crystal,flagFh=chi(args[0],h,k,l,T,Ex)
        result0=dtxrd1(thc,eta,phi,dc,Ex,P,crystal) # first iteration to find wh
        wh=result0[1]
        Ec=Eb*(1.0+wh)/sin(thc) # found center energy Ec
#             print "Ec = ", Ec
        crystal,flagFh=chi(args[0],h,k,l,T,Ec)
        result0=dtxrd1(thc,eta,phi,dc,Ec,P,crystal) # second iteration to find wh
        wh=result0[1]
        Ec=Eb*(1.0+wh)/sin(thc) # found center energy Ec
#             print "Ec = ", Ec
        # now determine constants at Ec, thc
        crystal,flagFh=chi(args[0],h,k,l,T,Ec)
        [[Chi0,Chih,Chih_],dh]=crystal
        [wh_s,wh,ep,dt,de,Tc,Rc]=dtxrd1(thc,eta,phi,dc,Ec,P,crystal)
        [eps_s,eps,eps_pr]=ep
        [dth_s,dth,dth_pr]=dt
        ################################################
        E1=-5.0*eps*Eb  # eV
        E2=5.0*eps*Eb   # eV
        Estep=(E2-E1)/1000.0 # eV
        ################################################
        Ex_v=Ec+arange(E1,E2,Estep)  # create energy vector

        Tplot=[]; Rplot=[]
        for Ex in Ex_v:
            result2=dtxrd1(thc,eta,phi,dc,Ex,P,crystal) # now calculate transmissivity/reflectivity
            Tplot=Tplot+[result2[5]]; Rplot=Rplot+[result2[6]]

        # convolution to calculate FWHM for the refl. curve
        if opts.conv!=0:
            fwhm0=0.1*eps*Ec # FWHM estimate for instrument resolution function
            #bw=0.5*fwhm0; b=[0.0,1.0,Ec,bw]; rf=lorentz(b,Ex_v)               # Lorentzian
            bw=0.5*fwhm0/sqrt(log(2.0)); b=[0.0,1.0,Ec,bw]; rf=gauss(b,Ex_v)   # Gaussian
            rf=rf/sum(rf)
            rc0=convolve(Rplot,rf,'same')
            stat0=curvestat(Ex_v,rc0,0.0)
            fwhm_ex=stat0[5]
        else:
            rc0=zeros(len(Rplot))
        #
        # plot all
        Eplot=1.0e3*(Ex_v-Ec) # meV
        outx0=1.0e-3*Ex_v     # keV
        outx1=Eplot           # meV
        outh0='E[keV] '       # header_col0
        outh1='E-Ec[meV] '    # header_col1

        plt.subplot(212)
        plt.plot(Eplot,Tplot,'b-')
        plt.ylabel('Transmissivity')
        plt.xlabel('E-E0 [meV]')
        #
        plt.subplot(211)
        plt.plot(Eplot,Rplot,'r-')
        if opts.conv!=0: plt.plot(Eplot,rc0,'k-')
        plt.ylabel('Reflectivity')

########
    elif flag=='e':
        Ex=1.0e3*float(args[9])
        Ec=Ex       # assume input energy as true central energy
        if Ex<Eb:
            fatalError('Ex must be greater than Eb')
        else:
            thb=arcsin(Eb/Ex)  # approximate angle from kinematic Braggs Law
            outFile.write('thb [deg] = '+str(r2d*thb)+' Bragg angle\n')

        if opts.polarization==1:
            P=cos(2.0*thb)
        crystal,flagFh=chi(args[0],h,k,l,T,Ec)
        result0=dtxrd1(thb,eta,phi,dc,Ec,P,crystal) # first iteration to find wh
        wh=result0[1]
        thc=arcsin((Eb/Ec)*(1.0+wh))        # determine central angle

        if opts.polarization==1:
            P=cos(2.0*thc)
        result0=dtxrd1(thc,eta,phi,dc,Ec,P,crystal) # second iteration to find wh
        wh=result0[1]
        thc=arcsin((Eb/Ec)*(1.0+wh))        # determine central angle

        if opts.polarization==1:
            P=cos(2.0*thc)
        # calculate constants at thc,Ec
        [[Chi0,Chih,Chih_],dh]=crystal
        [wh_s,wh,ep,dt,de,Tc,Rc]=dtxrd1(thc,eta,phi,dc,Ec,P,crystal)
        [eps_s,eps,eps_pr]=ep
        [dth_s,dth,dth_pr]=dt
        # thc - final iteration
        thc=arcsin((Eb/Ec)*(1.0+wh))
        #################################################
        th1=-5.0*dth  # rad
        th2=5.0*dth   # rad
        thstep= (th2-th1)/1000.0 # rad
        #################################################
        thb_v=thc+arange(th1,th2,thstep)  # create vector of angles

        Tplot=[]; Rplot=[]
        for thx in thb_v:
            if opts.polarization==1:
                P=cos(2.0*thx)
            result2=dtxrd1(thx,eta,phi,dc,Ec,P,crystal)  # now calculate transmissivity/reflectivity
            Tplot=Tplot+[result2[5]]; Rplot=Rplot+[result2[6]]

        # convolution to calculate FWHM
        if opts.conv!=0:
            fwhm0=0.1*dth  # FWHM estimate for instrument resolution function
            #bw=0.5*fwhm0; b=[0.0,1.0,thc,bw]; rf=lorentz(b,thb_v)          # Lorentzian
            bw=0.5*fwhm0/sqrt(log(2.0)); b=[0.0,1.0,thc,bw]; rf=gauss(b,thb_v)   # Gaussian
            rf=rf/sum(rf)
            rc0=convolve(Rplot,rf,'same')
            stat0=curvestat(thb_v,rc0,0.0)
            fwhm_th=stat0[5]
        else: rc0=zeros(len(Rplot))
        # plot results ##################################################
        thplot=1.0e6*(thb_v-thc) # urad
        outx0=thb_v*r2d          # deg
        outx1=thplot             # urad
        outh0='th[deg] '         # header_col0
        outh1='th-thc[urad] '    # header_col1

        plt.subplot(212)
        plt.plot(thplot,Tplot,'b-')
        plt.ylabel('Transmissivity')
        plt.xlabel('theta-theta0 [urad]')
        #plt.axis([1.0e6*th1,1.0e6*th2,0,auto])
        #
        plt.subplot(211)
        plt.plot(thplot,Rplot,'r-')
        if opts.conv!=0: plt.plot(thplot,rc0,'k-')
        plt.ylabel('Reflectivity')
        #plt.axis([1.0e6*th1,1.0e6*th2,0,2])

    else: fatalError('flag arg[8] either a (angle) or e (energy)')

###############################################################################
####### CONTINUE OUTPUT
###############################################################################

    G0=cos(thc)*sin(eta)*cos(phi)+sin(thb)*cos(eta)
    Gh=G0-2.0*sin(thb)*cos(eta)
    bh=G0/Gh
    #
    ## Calculate dispersion rate #################
    #Dr=2.0e3*sin(thc)*sin(eta)/(Ec*sin(thc-eta)) #murad/meV
    Dr=-2.0e3/(cos(phi)/tan(thc)-1.0/tan(eta))/Ec # now using (2.138)
    # notes:
    # this is angular dispersion rate in the dispersion plane
    # for phi=0 and phi=180 represents variation of dth_pr
    # it can't be maximized by choice of phi=90 since eta < thc
    # i.e. cot(eta) is always greater than cot(thc)

    outFile.write('bh = '+str(bh)+' asymmetry factor\n')
    outFile.write('##############################################################\n')
    outFile.write('# DYNAMICAL THEORY: \n')
    outFile.write('##############################################################\n')
    outFile.write('# Susceptibilities and Refraction Corrections: \n')
    outFile.write('chi_{0}  = ' +str(Chi0)+  '\n')
    outFile.write('chi_{h}  = ' +str(Chih)+  '\n')
    outFile.write('chi_{-h} = '+str(Chih_)+ '\n')
    outFile.write('wh(s) = '+str(wh_s)+  '\n')
    outFile.write('wh    = '+str(wh)+    '\n')
    outFile.write('#-------------------------------------------------------------\n')
    outFile.write('# Central Energy/Angle: \n')
    outFile.write('Ec [keV]  = '+str(1.0e-3*Ec)+'\n')
    outFile.write('thc [deg] = '+str(thc*r2d)+'\n')
    outFile.write('#-------------------------------------------------------------\n')
    outFile.write('# Energy and Angular Intriscis Widths (thick non-absorbing crystal):\n')
    outFile.write('eps_s  = ' +str(eps_s)+ '\n')
    outFile.write('eps    = ' +str(eps)+   '\n')
    outFile.write('eps_pr = ' +str(eps_pr)+'\n')
    outFile.write('Delta_E_s [meV]    = '+str(1.0e3*Ec*eps_s)+' \n')
    outFile.write('Delta_E [meV]    = '+str(1.0e3*Ec*eps)+' \n')
    outFile.write('Delta_E_pr [meV] = '+str(1.0e3*Ec*eps_pr)+' \n')
    outFile.write('dth_s [urad]  = '+str(1.0e6*dth_s)+ '\n')
    outFile.write('dth   [urad]  = '+str(1.0e6*dth)+   '\n')
    outFile.write('dth_pr [urad] = '+str(1.0e6*dth_pr)+ '\n')
    outFile.write('#-------------------------------------------------------------\n')
    outFile.write('# Additional Characteristics of the Reflection: \n')
    outFile.write('dE/dth [meV/urad] = '+str(1.0e-3*Ec/tan(thc))+' Tangent of the Braggs Law \n')
    outFile.write('Dr [urad/meV] = '+str(Dr)+' Intrinsic Angular Dispersion Rate \n')
    outFile.write('de [um] = '+str(1.0e-4*de)+' Extinction Length \n')
    #
    outFile.write('#-------------------------------------------------------------\n')

    if opts.conv!=0:
        outFile.write('# Widths Expected in Experiment (theory + inst.res.func.)\n')
        if flag=='a':
            outFile.write('Delta_E [meV] = '+str(1.0e3*fwhm_ex)+ '\n')
            outFile.write('inst. resolution [meV] = '+str(1.0e3*fwhm0)+ '\n')
        elif flag=='e':
            outFile.write('dth [urad] = '+str(1.0e6*fwhm_th)+ '\n')
            outFile.write('instrum. resolution [urad] = '+str(1.0e6*fwhm0)+ '\n')
        outFile.write('#-------------------------------------------------------------\n')
    #
    outFile.write('# Reflectivity/Transmissivity at Center: \n')
    outFile.write('Rc [%]  = '+str(Rc*100.0)+' Reflectivity \n')
    outFile.write('Tc [%]  = '+str(Tc*100.0)+' Transmissivity \n')
    outFile.close


    if opts.write is not None:
        header = '# '+'2beam '+str(args)+'\n' \
                 +'# version '+__version__+' by Stanislav Stoupin\n' \
                 +'# columns: \n' \
                 + '# '+outh0+outh1+'Tplot Rplot rc0\n'
        try:
            writeFile(opts.write,header,outx0,outx1,Tplot,Rplot,rc0)

        except IOError, e:
            fatalIOError(e)

    plt.show()


if __name__ == '__main__':
    main()
