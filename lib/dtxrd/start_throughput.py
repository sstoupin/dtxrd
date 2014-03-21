#!/usr/bin/env python

'''
a program to calculate throughput of a multicrystal configuration

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys
from scipy.special import legendre
from pylab import *

import os
if os.path.abspath(os.path.dirname(__file__)).split(os.sep)[-2] == 'lib':
    '''when running this script from the source directory'''
    sys.path.insert(0, os.path.abspath('..'))

from dtxrd.myio import *
from dtxrd.rotation import *
from dtxrd.curvestat import *
from dtxrd.chi import *
from dtxrd.thfind import *
from dtxrd.dtxrd2_k import *
from dtxrd.constants import *

__version__ = '0.16'
#-------------------------------------------------------------------------------------------
# SUMMARY OF CHANGES:
#-------------------------------------------------------------------------------------------
# 0.01 subroutine for iterative search of opt. angle
#      calculate throughput for monochromatic wave (all cryst. at opt. anlges)
# 0.02 rotation matrices + rules on assignment of axes
# 0.03 introduce transmission/reflection elements
# 0.04 implement search of thr - angle of exact backscattering
# 0.05 added pi polarization, added energy range as an argument, added integrated thruput
# 0.06 added rejection of Ex<Eb into 1 beam case (dtxrd_k v. 024)
# 0.07 integrated angular scans - can be called with -a option
# 0.08e bunch of tests to understand stuff - lots of thinking
# 0.09 angular scans work finally and understood how - see notes release 0.09
# 0.10 added option - last param crx[14] is now also angular shift if no scan is performed
# !!!  fixed angular condition for transition from a0 -> a0 when incident angle thr-th_in > 0.5*pi-thr
# 0.11 reintroduced entity card to preserve x-ray momentum during alignment !
# some angles should be equal e.g. in CDDW th5=th1_pr
# 0.12 averaged throughput and integrated ([meV]) throughput in the ouput
#      added "flat"(isotropic) angular distribution
# 0.13 can now take incoming photon energy distribution from file source_e.dat
# 0.14 introduce sapphire + new functional form (chi); now Ebragg -> Ebragg2
# 0.15 additional card ebb=9 to switch off alignment to backscattering crystal
#      added calculation of the cummulative dispersion rate
# 0.16 additional columns in the output file (Ex(keV) or th (deg))
#-------------------------------------------------------------------------------------------------------
#  things to do
# 1. vector dtxrd function?
# 2. define geometry by euler angles
# 3. introduce offsets in chi
# 4. take into account divergence in chi

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
'This program requires Optik, availible from http://optik.sourceforge.net/\n')

    USAGE = '%prog [OPTIONS...] func dpsi Ec dEx ne crystals_input_file \n'+'where:\n' \
    +'func - angular divergence distribution function for incoming - g(gaussian), l(lorentzian)\n' \
    +'dpsi - angular dispersion g - rms [urad] l - fwhm [urad] \n' \
    +'Ec - central energy [keV] \n' \
    +'dEx - energy range [meV] \n' \
    +'ne  - number of energy steps \n' \
    +'crystals_input_file - input file containing crystal parameters\n'
#
    VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>\n' \
    +'the program calculates throughput of a multicrystal system described in crystals_input_file'
    parser = OptionParser(usage=USAGE, version=VERSION)
#----------------------------------------------------------------------------------------------
    parser.add_option('-o', '--output', action='store', dest='output', default=None,
            help='write results to file F (defaults to stdout)', metavar='OUTFILE')
    parser.add_option('-a', '--angular_scan', type='int', action='store', dest='nth', default=0,
            help='do angular scan with nth points')
    parser.add_option('-s', '--source', type='int', action='store', dest='src', default=0,
            help='type of energy distribution for the source: 0 - flat (default), 1 - Cu K alpha, 9 - energy dist from file source_e.dat')
    parser.add_option('-w', '--write', action='store', dest='write', default=None,
            help='write data to file (default do not write)', metavar='DATAFILE')
    parser.add_option('-p', '--pi', action='store_const', const=1, dest='polarization', default=0,
            help='-p - pi polarization ')
#-----------------------------------------------------------------------------------------------
    opts, args = parser.parse_args(args)
    if len(args) != 6:
        parser.print_usage()
        sys.exit(1)

    try:
        args[1] = float(args[1])

    except ValueError:
        fatalError('dpsi must be a valid number')

    return opts, args


def main():
    opts, args = ParseArguments(sys.argv[1:])

    if opts.output is not None:
        try:
            outFile = open(opts.output, 'w')
        except IOError, e:
            fatalIOError(e)
    else:
        outFile = sys.stdout

######## read input parameters ######################################
    type=args[0]
    dpsi=1.0e-6*float(args[1])  # rms value - in my gauss sigma=sqrt(2)*dpsi then
    Ec=1.0e3*float(args[2]) # central energy converted to [eV]
    dEx=float(args[3])  # energy range to evaluate [meV]
    ne=int(args[4]) # number of steps in energy
    nth=opts.nth
######## IF FIXED ENERGY SOURCE - ASSIGN ENERGY #######################################################
    if opts.src==1:   # Cu Kalpha1
        Ec=8047.84
        print "performing angular alignment at Cu Kalpha1 energy [keV] ", Ec*1.0e-3
    #
######## READ CRYSTAL DATA ############################################################################
    filein=args[5]
    header,crdata=readCrystal(filein)
    items=range(0,len(crdata),1)
    pairs=range(0,len(crdata)-1,1)
    # check if backscattering card is on (equal to 1) for any crystal
    # and if true - perform angle angular tuning at the backscattering energy of the first crystal
######### for which the ebb card is 1 #################################################################
    ebbcard=0
    for crx in crdata:
        if int(crx[11])== (1 or 2) and ebbcard==0:  # if the first EBB condition is found
            el=crx[0]
            h=float(crx[1])
            k=float(crx[2])
            l=float(crx[3])
            eta=float(crx[4])/r2d
            phi=float(crx[5])/r2d
            T=float(crx[6])
            dc=1.0e7*float(crx[7])
#              ddth=float(crx[12])
            # do stuff to find true Er
            thr,Er = thr_find(el,h,k,l,eta,phi,T,dc,Ec)
            #if opts.polarization==1:
            #          thb=arcsin(Eb/Ec)
            #          P=cos(2.0*thb)
            #else: P=1.0
            #-------------------------------
            print "thr = ", thr*r2d
            print "Er [keV] = ", Er*1.0e-3
            Ec=Er
#              th_in=thr-1.0e-6*ddth   ; print "th_in = ", th_in*r2d#
#              if opts.polarization==1: P=cos(2.0*th_in)
#              result0=dtxrd(el,h,k,l,th_in,eta,phi,a,dh,T,dc,Ec,P)
#              wh=result0[4]
#              Ec=Eb*(1.0+wh)/sin(th_in)
            print "performing angular optimization at energy Ec [keV] = ", Ec*1.0e-3
            ebbcard=1
####################################################################################
######## FILL LISTS TO STORE CRYSTAL PARAMETERS
####################################################################################
    el_l=[]; h_l=[]; k_l=[]; l_l=[]; eta_l=[]; phi_l=[]; T_l=[]; dc_l=[]; ss_l=[]
    crystal=[]; thc_l=[]; thc_pr_l=[]; eps_l=[]; bh_l=[]
    Tpt0=1.0; Ty=[]; Ry=[] #;Chi_=[]
    flag_l=[]; ddth_l=[]; ebb_l=[] ;entity_l=[]
    cntas=0; cntbb=[]
    Dr_l=[]
    #ans_l=[]
    if ebbcard==0:
        print "performing angular optimization at energy Ec [keV] = ", Ec*1.0e-3
    cr_cnt=0; cr_cnt_l=[cr_cnt]
    for crx in crdata:
        el=crx[0];                el_l=el_l+[el]
        h=float(crx[1]);          h_l=h_l+[h]
        k=float(crx[2]);          k_l=k_l+[k]
        l=float(crx[3]);          l_l=l_l+[l]
        eta=float(crx[4])/r2d;    eta_l=eta_l+[eta]   #; print "cos(eta) =", cos(eta)
        phi=float(crx[5])/r2d;    phi_l=phi_l+[phi]
        T=float(crx[6]);          T_l=T_l+[T]
        dc=1.0e7*float(crx[7]);   dc_l=dc_l+[dc]      # [A] crystal thickness
        ss=int(crx[8]);           ss_l=ss_l+[ss]      # sign of crystal configuration (multiplier)
        #------------------------------------------------------------------------------------------
        flag=crx[9];
        if flag!='R' and flag!='T':
            fatalIOError("check crystal configuration - must be either R(reflection) or T(transmisson)")
        else: flag_l=flag_l+[flag] # R - reflection T - transmission
        #
        entity=int(crx[10]);           entity_l=entity_l+[entity] # something that exists as a unit
        ebb=int(crx[11]);              ebb_l=ebb_l+[ebb]
        ddth=float(crx[12]);           ddth_l=ddth_l+[ddth]
        ans=int(crx[13])
        anr=float(crx[14])
#           if nth!=0:
#              ans=int(crx[13])            #;ans_l=ans_l+[ans]
#              if ans!=0:
#                 anr=float(crx[14])
        if nth==0:
            ans=0
            #anr=0.0
        print "Crystal "+str(entity)+' , structure '+str(el)
        ###################################################################################################################
        ##  Find optimal angle at the given energy
        ###################################################################################################################
        crystalx,flagFh=chi(el,h,k,l,T,Ec); crystal=crystal+[crystalx]
        if flagFh==0:
            fatalError('forbidden reflection: structure amplitude |exp(iHr)| for the chosen set of Miller indicies is too small (< 1e-6)')
#           Chixx=crystalx[0]
        # rangle of angles for angular scan
        if ans!=0:
            cntas=cntas+1
            if ans==2:
                hnth=int(0.5*float(nth))
                hr=logspace(0.01,1,hnth,base=100)-1.0
                hr_pos=append(0.0,hr)
                hr_neg=sort(-hr)
                thva=1.0e-6*append(hr_neg,hr_pos)*float(hnth)*anr/99.0
            else:
                thva=1.0e-6*arange(-0.5*float(nth)*anr,0.5*float(nth)*anr+anr,anr) # create vector for angular scan
            #
            if cntas==1 and opts.nth !=0:
                thscan=thva
                cr_cnt_ans=cr_cnt
        else:
            thva=1.0e-6*anr*ones(nth+1)
        #
        if opts.polarization==1:
            dh=crystalx[1]; Eb=0.5*hpl*cl/dh
            thb=arcsin(Eb/Ec)
            P=cos(2.0*thb)
        else: P=1.0
#           Pv=1.0+zeros(nth+1)
#!!!
        if ebb==0:
#              thcx,thcx_pr,epsx,bhx,Tyx,Ryx=thmax_find(Ec,eta,phi,dc,crystalx,P,'R')
            thcx,thcx_pr,epsx,bhx,Tyx,Ryx=thc_find(Ec,eta,phi,dc,crystalx,P)
            th_in=thcx-1.0e-6*ddth
            Dr=-2.0e3*ss/(cos(phi)/tan(th_in)-1.0/tan(eta))/Ec # disp. rate using (2.138)
            ####################################################################
            ## ENTITY ##########################################################
            ####################################################################
#              if cr_cnt>1:
#               tmp=[]
#               for (x,i) in map(None,entity_l,cr_cnt_l):
#                if tmp.__contains__(x):
#                  gamma=pi-mean(thc_pr_l[i-1])-mean(thc_l[i-1])
#                  print 'gamma [urad] = ', gamma*1.0e6
#                  th_in=mean(thc_pr_l[i-2])-float(ss_l[i-1])*gamma
#                else: tmp.append(x)
            ####################################################################
            ## CALCULATE th_out ################################################
            ####################################################################
            if opts.polarization==1: P=cos(2.0*th_in)
            k_in=[cos(th_in)*cos(phi),cos(th_in)*sin(phi),-sin(th_in)]
            k_out,thcxx,dthxx,epsxx,bhxx,Tyxx,Ryxx=dtxrd0_k(k_in,eta,dc,Ec,P,crystalx)
            thcxx_pr=arcsin(k_out[2])
            # vectorize ------------------------------------------------------------------------------
            thv_in=th_in-thva
            thv_out=thcxx_pr+thva    # because it is rotation of crystal plane not change angle or ray
                                     # and the reflected angle is from opposite side
        #----------------------------------------------------------------------------
        elif ebb==1 or ebb==2 or ebb==9:
            thr,Er=thr_find(el,h,k,l,eta,phi,T,dc,Ec)
            thcx=thr;                   print "thr = ", thr*r2d
            bthr=0.5*pi-thr;            print "bthr [urad] = ", bthr*1.0e6
            th_in=thcx-1.0e-6*ddth
            Dr=-2.0e3*ss/(cos(phi)/tan(th_in)-1.0/tan(eta))/Ec # disp. rate using (2.138)
            if opts.polarization==1: P=cos(2.0*th_in)
            k_in=[cos(th_in)*cos(phi),cos(th_in)*sin(phi),-sin(th_in)]
            k_out,thcxx,dthxx,epsxx,bhxx,Tyxx,Ryxx=dtxrd0_k(k_in,eta,dc,Ec,P,crystalx)
            thcxx_pr=arcsin(k_out[2]) #; print "thcxx_pr = ", thcxx_pr*r2d
            #--------------------------------
            # Reflected big theta tilda (from virtual planes)
            oth_red=(pi-thr-thcxx_pr); print "oth_red [urad] ", 1.0e6*oth_red
            oth_blue=(thr-thcxx_pr); print "oth_blue [urad] ", 1.0e6*oth_blue
            #--------------------------------
            # Choosing whether it is beyond bthr or not
            # only need to do this for the central angle which is dtxrd calculated!!!
            bth_t=abs(thr-th_in)
            if bth_t < bthr:
                th_out=pi-thcxx_pr
            else:
                th_out=thcxx_pr
            # vectorize -----------------------------------------------------------------------------
            thv_in=th_in-thva
            thv_out=th_out+thva # because it is rotation of crystal plane not change angle of ray
                                # and the reflected angle is from opposite side
            #----------------------------------------------------------------------------------------
        else:
            fatalIOError("check crystal configuration file - ebb either 0,1,2 or 9 (no alignment to ebb)")
##########################################################################################################
########## FILL LISTS (VECTORS CORRESPONDING TO DIFFERENT ANGLES OF THE ROTATING CRYSTAL(S) ##############
##########################################################################################################
        if flag=='T':
            thv_out=thv_in
            bhxx=1.0
            Dr=0.0
        # fill the lists
        thc_l=thc_l+[thv_in] #; print thc_l
        thc_pr_l=thc_pr_l+[thv_out] #; print thc_pr_l
        eps_l=eps_l+[epsxx]; bh_l=bh_l+[bhxx]
        Ty=Ty+[Tyxx]; Ry=Ry+[Ryxx]
        Dr_l=Dr_l+[Dr]  # dispersion rates
#           Chi_=Chi_+[Chixx]
        print "thva [deg] ", thva*r2d
        cr_cnt=cr_cnt+1; cr_cnt_l=cr_cnt_l+[cr_cnt]
#################################################################################
########### ENTITY ##############################################################
#################################################################################
    tmp=[]
    for (x,i) in map(None,entity_l,items):
        if tmp.__contains__(x):
            gamma=pi-mean(thc_pr_l[i-1])-mean(thc_l[i-1])
            if flag_l[i-2]=='T': gamma=-gamma; print "YES"   # to account for T->ebb->R case
            print 'gamma [urad] = ', gamma*1.0e6
            th_in=mean(thc_pr_l[i-2])-float(ss_l[i-1])*gamma
            if opts.polarization==1: P=cos(2.0*thv_in)
            k_in=[cos(th_in)*cos(phi_l[i]),cos(th_in)*sin(phi_l[i]),-sin(th_in)]
            result=dtxrd0_k(k_in,eta_l[i],dc_l[i],Ec,P,crystal[i])
            k_out=result[0]
            th_out=arcsin(k_out[2])
            thva=mean(thc_l[i])-thc_l[i]
            thc_l[i]=th_in-thva
            if flag_l[i]=='T': thc_pr_l[i]=thc_l[i]
            else:              thc_pr_l[i]=th_out+thva
        else: tmp.append(x)
###############################################################################################
########  CALCULATE THROUGHPUT IN OPTIMAL GEOMETRY ############################################
###############################################################################################
    for i in items:
        if flag_l[i]=='T': Tpt0=Tpt0*Ty[i]; print "Ty = ", Ty[i]
        else:              Tpt0=Tpt0*Ry[i]; print "Ry = ", Ry[i]

####### find total energy acceptance of the first crystal #####################################
    th1=mean(thc_l[0]); eta1=eta_l[0]; phi1=phi_l[0]; dc1=dc_l[0]
    eps1=eps_l[0]; dE1=1.0e3*Ec*eps1; print "Energy acceptance region of cr1 [meV] = ", dE1
    dedth=1.0e-3*Ec/tan(th1) #meV/urad
    dE1psi=dedth*1.0e6*dpsi  # meV
    dE1tot=sqrt(dE1**2.0+dE1psi**2.0); print "Energy acceptance total cr1 [meV] = ", dE1tot
###############################################################################################
## create array of energies and angular divergence
###############################################################################################
    Eran=dEx # [meV]
#        step=2.0/ne
#        Ev=Ec+Eran*1.0e-3*arange(-1.0,1.0+step,step)  # [eV] for calculation purposes
    pe=legendre(ne); pe_w=transpose(array(pe.weights))
    Ev=Ec+1.0e-3*Eran*pe_w[0]
    Ew=pe_w[1]
#
    if opts.src==0:
        edist=ones(len(Ev))
        Norm_E=1.0
    elif opts.src==1:
        from beams import CuKa, CuKa_hartwig_sc
        edist=CuKa_hartwig_sc(Ev)
    elif opts.src==9:
        from scipy.interpolate import interp1d
        d1,d2=readFile('source_e.dat')
        eev=d2[:,0]
        rrv=d2[:,1]
        rrf=interp1d(eev,rrv)
        edist=rrf(Eran*pe_w[0])
#           figure(99)
#           plot(eev,rrv,'k-')
#           plot(Eran*pe_w[0],edist,'b-')
#           show()
#
    Norm_E=0.5*sum(Ew*edist); print "Norm_E = ", Norm_E
####### angular divergence ###########################################################################################
    if type !='1r':
        npsi=200
        pa=legendre(npsi); pa_w=transpose(array(pa.weights))
        thv=3.0*dpsi*pa_w[0]  # 3.0*dpsi=0.5*(b-a) integration from -3*dpsi to 3*dpsi; for Gauss +-3*sigma - 99.7%
        thw=pa_w[1]
    #thva=1.0e-6*arange(-0.5*nth*ddth,0.5*nth*ddth+ddth,ddth)  # create vector of incident angles for angular scan
    if type =='l':
        bw=0.5*dpsi # lorentz
        param_l=[0.0,1.0,th1,bw]
        psiv=lorentz(param_l,th1+thv)
        fwhm_psi=dpsi
    elif type == 'g':
        sigma=sqrt(2.0)*dpsi #; print sigma
        param_g=[0.0,1.0,th1,sigma]
        psiv=gauss(param_g,th1+thv)
        fwhm_psi=2.0*sigma*sqrt(log(2.0))
    elif type =='1r':
        npsi=1
        psiv=ones(1)
        thv=array([0.0])
        fwhm_psi=0.0
        thw=2.0*ones(1)
    elif type =='f':     #flat from -3*dpsi to 3*dpsi
        fwhm_psi=float('inf')
        psiv=ones(len(thv))
    else:
        fatalError('res. function type is either "l" ,"g", "1r" or "f" at the moment')

    Norm_th=0.5*sum(thw*psiv); print "Norm_th = ", Norm_th

############################################################################################################
## Throughput calculation
############################################################################################################
    # calculate transformation matrix
    tr=[]
    ass=[]
    for ssx in ss_l:
        assx=ones(nth+1)*ssx
        ass=ass+[assx]

    thc_a=array(thc_l)*array(ass) # array of thc values taken with an appropriate crystal sign
    thc_pr_a=array(thc_pr_l)*array(ass)
    # DIAGNOSTICS
    for (th_in,th_out) in map(None,thc_a,thc_pr_a):
        print th_in*r2d, th_out*r2d

#        print thc_a*r2d
#        print thc_pr_a*r2d
    ####################################################################################################
    # find direction of y-axis
    ####################################################################################################
#        ydir=ss_l
    ydir=[0]*len(crystal) # positive inboard, corresponds to + counterclockwise rotation
    for i in items:
        if   abs(phi_l[i]) < 0.5*pi: # aka == 0  calculating phi dependencies will require rot_z by that phi
            if ss_l[i]==1:
                ydir[i]=1
            else:          ydir[i]=-1
        elif abs(phi_l[i]) > 0.5*pi: # aka == 180
            if ss_l[i]==1: ydir[i]=-1
            else:          ydir[i]=1
        else:              ydir[i]=1
    print ydir

    for i in pairs:
        #---------------------------------------------------
        if flag_l[i]=='T':
#              dthxx=-thc_a[i]+pi-thc_a[i+1]
            dthxx=(-thc_a[i]+thc_a[i+1])
        else:
            dthxx=(thc_pr_a[i]+thc_a[i+1]) #; print dthxx*r2d

        if ebb_l[i+1]==2:
            if ss_l[i+1]==-1:
                dthxx=-(thc_pr_a[i]+pi-thc_a[i+1])
        #---------------------------------------------------
        if ss_l[i]*ss_l[i+1]==-1:
            dthxx=dthxx+pi
            if ydir[i]*ydir[i+1]==-1: bonus=rot_z(pi)
            else:                     bonus=identity3
        else:
            if ydir[i]*ydir[i+1]==-1: bonus=inv_x
            else:                     bonus=identity3

        tr_xx=[]
        for t in dthxx:
            if ydir[i]==1:
                try_xx=rot_y(t)
            else:
                try_xx=rot_y(-t)
            tr_xx=tr_xx+[dot(bonus,try_xx)]
#            print tr_xx
        tr=tr+[tr_xx]
    print "------------------------------------------------"
################################################################################################################
###  MAIN LOOP
################################################################################################################
    if opts.nth == 0:
        thscan=[0]
    if opts.write is not None:
        if opts.nth !=0:
            header = '# '+'throughput '+str(args)+'\n' \
                     +'# version '+__version__+' by Stanislav Stoupin <sstoupin@gmail.com>\n' \
                     +'# columns: (th-thc)[urad] Tptt th[deg]\n' \
                     + '# \n'
            fig0 = plt.figure(0)
            ion()
            plt.ylabel('Througput [norm.u.]')
            plt.xlabel('angle [urad]')

        else:
            header = '# '+'throughput '+str(args)+'\n' \
                     +'# version '+__version__+' by Stanislav Stoupin <sstoupin@gmail.com>\n' \
                     +'# columns: (E-Ec)[meV] Tpte Ex[keV]\n' \
                     + '# \n'
            plt.xlabel('Energy [meV]')
#-----------------------------------------------------------------------------------------------
    tra=[[row[i] for row in tr] for i in range(nth+1)]  # transpose list of matrixes tr
#
    Tptt=[]
    thscan0=[]
#        count=0
    for (thax,thscanx,trax) in map(None,thc_l[0],thscan,tra):   # ANGULAR SCAN
#         print "thax = ", thax*r2d
#         print "trax " , trax
        print "rocking angle [urad] ", thscanx*1.0e6
        Tpte=[]
        Eplot=[]
        for (ex,edx) in map(None,Ev,edist):
        #-----------------------------------------------------------------------------------
            tpty=[]
            for (thx,psix) in map(None,thv,psiv): # take each ray and propagate
#                print "thx = ", thx
                k1=[cos(thx+thax)*cos(phi1),cos(thx+thax)*sin(phi1),-sin(thx+thax)]
                if opts.polarization==1:
                    P=cos(2.0*(thx+thax))
                result1=dtxrd1_k(k1,eta1,dc1,ex,P,crystal[0])
#                result1=dtxrd_k(k1,eta1,dc1,ex,P,crystal[0])
                if flag_l[0]=='T': Tptx=result1[1]*psix*edx; kout=k1
                else:              Tptx=result1[2]*psix*edx; kout=result1[0]
#                print "cryst "+str(1)+" E: "+str(1.0e3*(ex-Ec))+" out. angle ", arcsin(kout[2])*r2d
                #---------------------------------------------------------------
                for i in pairs:
                    kin=dot(trax[i],kout)
                    #DIAGNOSTICS ###########################################################################
                    #if i==3:
                    #  print "pair: "+str(i)+" E: "+str(1.0e3*(ex-Ec))+" inc. angle ", -arcsin(kin[2])*r2d
                    #if i==4:
                    #  print "pair: "+str(i)+" E: "+str(1.0e3*(ex-Ec))+" inc. angle ", -arcsin(kin[2])*r2d
                #########################################################################################
                    if opts.polarization==1:
                        P=cos(2.0*(arcsin(kin[2])))
                    result=dtxrd1_k(kin,eta_l[i+1],dc_l[i+1],ex,P,crystal[i+1])
                    #result=dtxrd_k(kin,eta_l[i+1],dc_l[i+1],ex,P,crystal[i+1])
                    if flag_l[i+1]=='T': Tptx=Tptx*result[1]; kout=kin
                    else:                Tptx=Tptx*result[2]; kout=result[0]
                    ########################################################################################
#                    if i==0:
#                      print "pair: "+str(i)+" E: "+str(ex-Ec)+" out. angle ", arcsin(kout[2])*r2d
                #---------------------------------------------------------------------------------
                tpty=tpty+[Tptx]
            #-------------------------------------------------------------------------------------
            Tpte=Tpte+[0.5*sum(thw*tpty)/Norm_th]  # normalized by interval 2*thvmax
            #if opts.nth == 0:
            Eplot=Eplot+[1.0e3*(ex-Ec)]
#
        Tpt_int=0.5*sum(Ew*Tpte)/Norm_E           # normalized by interval 2*Eran
        Tptt=Tptt+[Tpt_int]
        thscan0=thscan0+[thscanx]
#         count=count+1

#############################################################################################
## OUTPUT
#############################################################################################
        if opts.write is not None:
            if opts.nth !=0:
                try:
                    thplot=(mean(thc_l[cr_cnt_ans])+array(thscan0))*r2d
                    writeFile(opts.write,header,1.0e6*array(thscan0),array(Tptt),thplot)
                except IOError, e:
                    fatalIOError(e)
                #-------------------------------------------
                # progress: plot angular scan at each angle
                #-------------------------------------------
                fig0.clear()
                plt.plot(1.0e6*array(thscan0),array(Tptt),'b-')
                #-------------------------------------------
                # diagnostics: plot thruput at each angle
                #-------------------------------------------
                #plt.xlabel('Energy [meV]')
                #plt.semilogy(Eplot,Tpte, label=str(thscanx))
                #plt.axis([-Eran,Eran,1e-20,1.0])
                #plt.legend(loc='upper right')
                #-------------------------------------------
                plt.draw()
            else:
                try:
                    writeFile(opts.write,header,array(Eplot),array(Tpte),1.0e-3*array(Ev))
                except IOError, e:
                    fatalIOError(e)
                #
#               plt.plot(array(Eplot),array(Tpte),'b-')
#               plt.draw()
    ##########################################################################
    ## Calculate cummulative dispersion rate
    ##########################################################################
    Din=0.0
    print "asym. factors: ", bh_l
    print "Dispersion rates [urad/meV]: ", Dr_l
    for (bh,Dr) in map(None,bh_l,Dr_l):
        Dout=bh*Din+Dr
        print "Dout = ", Dout
        Din=Dout

######  BEGIN OUTPUT #############################################################
    outFile.write('##############################################################\n')
    outFile.write('##### THRUPUT v'+__version__+' ########################################\n')
    outFile.write('##### Author: Stanislav Stoupin ## sstoupin@aps.anl.gov ######\n')
    outFile.write('##############################################################\n')
    outFile.write('Throughput (monochromatic ray) = '+str(Tpt0)+'\n')
    outFile.write('Incident divergence [urad] rms = '+str(1.0e6*dpsi)+'\n')
    outFile.write('Incident divergence [urad] fwhm = '+str(1.0e6*fwhm_psi)+'\n')
    outFile.write('Cummulative dispersion rate [urad/meV] = '+str(Dout)+'\n')
#        outFile.write('Throughput (divergent ray) = '+str(Tpt_dpsi)+'\n')
    outFile.write('Averaged throughput  = '+str(Tpt_int)+'\n')
    outFile.write('Integrated throughput [meV] = '+str(2.0*Eran*Tpt_int)+'\n')
#        outFile.write('--------------------------------------------------------------\n')
######################################################################################################

#        if opts.write is not None:
#               try:
#                    writeFile(opts.write,header,array(Eplot),array(Tpte))
#               try:
#                    writeFile(opts.write,header,1.0e6*array(thscan),array(Tptt))
#               except IOError, e:
#                    fatalIOError(e)
#######################################################################################################
# Plot the result
#######################################################################################################
####### log Y ###################################
    f1=plt.figure(101)
    if opts.nth !=0:
        plt.semilogy(1.0e6*thscan,Tptt,'r-')
    else:
        plt.semilogy(Eplot,Tpte,'r-')
        plt.axis([-Eran,Eran,1e-5,1.0])
######## linear Y ##############################
    f2=plt.figure(102)
    if opts.nth !=0:
        plt.plot(1.0e6*thscan,Tptt,'r-')
    else:
        plt.plot(Eplot,Tpte,'r-')
#
    plt.show()

#######  CURVE STATISTICS ##############################################################################
#        try:
    if opts.nth !=0:
        stat1=curvestat(array(1.0e6*thscan),array(Tptt),0.0)
        rmax=stat1[1]
        fwhm_t=stat1[5]
        outFile.write('FWHM [urad] = '+str(fwhm_t)+'\n')
        outFile.write('Rmax = ' +str(rmax)+'\n')
    else:
        stat1=curvestat(array(Eplot),array(Tpte),0.0)
        E0=1.0e-3*(Ec+1.0e-3*stat1[0])
        rmax=stat1[1]
        fwhm_e=stat1[5]
        outFile.write('FWHM [meV] = '+str(fwhm_e)+'\n')
        outFile.write('Emax [keV] = '+str(E0)+'\n')
        outFile.write('Rmax = '+str(rmax)+'\n')
#        except ValueError:
#           fatalError('FWHM not found - failed to find peak')

    outFile.close


if __name__ == '__main__':
    main()
