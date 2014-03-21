#!/usr/bin/env python

'''
a program to plot and calculate parameters of a reflectivity curve

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys
#from numpy import *
from pylab import *
#from scipy import *
from scipy.optimize import *
from scipy.interpolate import interp1d

from myio import *
from curvestat import *
from deriv import *

__version__='0.27'
proginfo = 'peak v' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>'

########################################################################
# modifications:
########################################################################
# v0.2
#  1. implemented background subtraction before estimation of fwhm
#  2. implemented linear interpolation for finding th_neg and th_pos
#  3. combined methods into subroutine curvestat
#  4. optional: user can define background value (default: min(r))
# v0.21
#  1. got rid of xpd (myio instead) and
#  2. got rid of Scientific.Func..python-scientific
#  3. added common plot
# v0.22
#  1. included statistics
# v0.23
#  1.various cosmetic changes
# v0.24
#  1. add derivative calculation and plotting
# v0.25 add multiply option (e.g. invert data)
# v0.26 introduce contrast as an option
# v0.27 add -n option to normalize for common plot
########################################################################
# Todo list
# 1. interpolate data
# 2. background function
#---------------------------------------------------
# Functions
#---------------------------------------------------

cbank=('b','g','r','c','m','y')

def residuals1(p,y,x):
    return  y-gauss(p,x)

def residuals2(p,y,x):
    return y-lorentz(p,x)
#############################################################################################

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

    USAGE = '%prog [OPTIONS...] FILENAME1 FILENAME2'
    VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>\n' \
    +'the program performs characterization of reflectivity curve\n' \
    +'returning various peak and curve width parameters'
    parser = OptionParser(usage=USAGE, version=VERSION)
#################################################################################################
    # Behavior
#################################################################################################

    #Output
    parser.add_option('-o', '--output',
            action='store',
            dest='output',
            default=None,
            help='write results to file F (defaults to stdout)',
            metavar='F')
    parser.add_option('-p', '--plot',
            action='store_const',
            const=1,
            dest='plotting',
            default=0,
            help='plot results')
    parser.add_option('-c', '--contrast',
            action='store_const',
            const=1,
            dest='cntr',
            default=0,
            help='calculate spectral contrast at 1.578 fwhm')
    parser.add_option('-u', '--unit',
            action='store',
            type='float',
            dest='unit',
            default='1',
            help='conversion factor for x-axis (th)',
            metavar='unit')
    parser.add_option('-m', '--mult',
            action='store',
            type='float',
            dest='mult',
            default='1.0',
            help='conversion factor for y-axis (r)',
            metavar='mult')
    parser.add_option('-x', '--name_x',
            action='store',
            dest='nx',
            default='t [raw u.]',
            help='x-axis label/unit',
            metavar='F')
    parser.add_option('-y', '--name_y',
            action='store',
            dest='ny',
            default='r [raw u.]',
            help='y-axis label/unit',
            metavar='F')
    parser.add_option('-b', '--bkg',
            action='store',
            type='float',
            dest='bkg',
            default='-1',
            help='user defined background for calculation of fwhm',
            metavar='bkg')
    parser.add_option('-t','--th',
            action='store',
            type='string',
            dest='n_th',
            default='1',
            help='th values are in column number n_th',
            metavar='n_th')
    parser.add_option('-r','--r',
            action='store',
            type='string',
            dest='n_r',
            default='2',
            help='intensity values are in column number n_r',
            metavar='n_r')
    parser.add_option('-d', '--deriv',
            action='store',
            type='int',
            dest='nder',
            default='-1',
            help='interpolation + derivative',
            metavar='nder')
    parser.add_option('-n', '--normalize',
            action='store_const',
            const=1,
            dest='norm',
            default=0,
            help='normalize for plotting',
            metavar='norm')
    ##########################################################################################
    #Parameters
    ##########################################################################################
    opts, args = parser.parse_args(args)
    if len(args) < 1:
        parser.print_usage()
        sys.exit(1)
#        elif len(args) > 1:
#                parser.print_usage()
#                sys.exit(1)
    return opts, args
##################################################################################################
## Main stuff
##################################################################################################


def main():
    opts, args = ParseArguments(sys.argv[1:])

    if opts.output is not None:
        try:
            outFile = open(opts.output, 'w')
        except IOError, e:
            fatalIOError(e)
    else:
        outFile = sys.stdout
    #--------------------------------------------------------
    # initialize variables :
    #--------------------------------------------------------

    count=0
    cf=opts.unit
    nder=int(opts.nder)
    outFile.write('conversion factor for th units = '+str(cf)+'\n')
    #--------------------------------------------------------
    # do stuff
    #--------------------------------------------------------
    fwhm_=[]; com_=[]
#        c_pos=[]; c_neg=[]

    for fn in args:

        outFile.write('filename: '+fn+'\n')
        count=count+1
        d1,d2=readFile(fn)
        th=d2[:,int(opts.n_th)-1]
        r=d2[:,int(opts.n_r)-1]
        N1=len(th); N2=len(r)
        mult=float(opts.mult)
        r=mult*r
        # interpolation ########################
        rint=interp1d(th,r)
        #
        nint=1000
        th1=th[0]
        th2=th[len(th)-1]
        step=abs(th2-th1)/float(nint)
        thint=arange(th1,th2+step,step)
        ############################################
        # derivative
        ############################################
        if nder > 0:
            th_orig=th
            r_orig=r
            th,r=derivs(th,r,nder)
#----------------------------------------------------------------
#            k1=1
#            k2=len(temp1)
#            temp2=temp2[k1:k2]
#            alpha=alpha[k1:k2]
#            utexp=utexp[k1:k2]
#            DT=DT[k1:k2]
#----------------------------------------------------------------

        bkg=0.25*(r[0]+r[1]+r[len(r)-2]+r[len(r)-1])
        result=curvestat(th,r,bkg)
        th_max=result[0]; r_max=result[1]
        th_neg=result[2]; th_pos=result[3]; th_mid=result[4]
        fwhm=result[5]; com=result[6]
        if opts.cntr !=0:
        # calculate contrast
            thc_pos=th_mid+1.578*fwhm
            thc_neg=th_mid-1.578*fwhm
            c_pos=rint(th_mid)/rint(thc_pos)
            c_neg=rint(th_mid)/rint(thc_neg)

        fwhm_=fwhm_+[fwhm]; com_=com_+[com]
#            c_pos_=c_pos_+[c_pos]; c_neg_=c_neg_+[c_neg]

        # reference in theta
        if fn==args[0]:
            th0_mid=th_mid
  #--------------------------------
  # gaussian fit
  #--------------------------------
        p_init=[bkg, (r_max-bkg), th_mid, fwhm/sqrt(4.0*log(2.0))]
        x=th; y=r
        gaussfit=leastsq(residuals1, p_init, args=(y,x), full_output=1)
        p_fin=gaussfit[0]
        r_g=gauss(p_fin,thint)
        bkg_g=p_fin[0]; r_max_g=p_fin[1]; th_max_g=p_fin[2]; fwhm_g=sqrt(4.0*log(2.0))*p_fin[3]
#            fwhmg_=fwhmg_+[fwhm_g]
#            thmaxg_=thmaxg_+[th_max_g]

        chisq = abs(sum((residuals1(p_fin,y,x))**2.0))
        r2 = sum(residuals1(p_fin,y,x)**2.0)/sum((r_g-p_fin[0])**2.0)
        stdv = sqrt(1.0/N1*chisq)
        corr = stdv*sqrt(abs(gaussfit[1]))
        uncert = diag(corr)

        outFile.write('------------------------------------------------------------\n')
        outFile.write('Gaussian fit: \n')
        outFile.write('fwhm_g = '+str(cf*fwhm_g)+'\n')
        outFile.write('th_max_g = '+str(cf*th_max_g)+'\n')
        outFile.write('r_max_g = '+str(r_max_g)+'\n')
        outFile.write('bkg_g = '+str(bkg_g)+'\n')
        outFile.write('chisq = '+str(chisq)+'\n')

  #--------------------------------
  # lorentzian fit
  #--------------------------------
        p_init=[bkg, (r_max-bkg), th_mid, 0.5*fwhm]
        x=th; y=r
        lorentzfit=leastsq(residuals2, p_init, args=(y,x), full_output=1)
        p_fin=lorentzfit[0]
        r_l=lorentz(p_fin,thint)
        bkg_l=p_fin[0]; r_max_l=p_fin[1]; th_max_l=p_fin[2]; fwhm_l=2.0*p_fin[3]  #correct! p3=1/2*Gamma

        chisq = abs(sum((residuals2(p_fin,y,x))**2.0))
        r2 = sum(residuals1(p_fin,y,x)**2.0)/sum((r_l-p_fin[0])**2.0)
        stdv = sqrt(1.0/N1*chisq)
        corr = stdv*sqrt(abs(gaussfit[1]))
        uncert = diag(corr)

        outFile.write('-------------------------------------------------------------------\n')
        outFile.write('Lorentzian fit: \n')
        outFile.write('fwhm_l = '+str(cf*fwhm_l)+'\n')
        outFile.write('th_max_l = '+str(cf*th_max_l)+'\n')
        outFile.write('r_max_l = '+str(r_max_l)+'\n')
        outFile.write('bkg_l = '+str(bkg_l)+'\n')
        outFile.write('chisq = '+str(chisq)+'\n')

     #----------------------------------
     # refine fwhm - optional
     #----------------------------------
        if opts.bkg !=-1:
            result=curvestat(th,r,float(opts.bkg))
            th_max=result[0]; r_max=result[1]
            th_neg=result[2]; th_pos=result[3]; th_mid=result[4]
            fwhm=result[5]; com=result[6]

  #--------------------------------
  # PRINT RESULTS
  #--------------------------------
        outFile.write('------------------------------------------------------------\n')
        outFile.write('Generic FWHM\COM: \n')
        outFile.write('fwhm = '+str(cf*fwhm)+'\n')
        outFile.write('th_mid = '+str(cf*th_mid)+'\n')
        outFile.write('th_max = '+str(cf*th_max)+'\n')
        outFile.write('r_max = '+str(max(r))+'\n')
        outFile.write('COM = '+str(cf*com)+'\n')
        outFile.write('bkg = '+str(bkg)+'\n')
        if opts.cntr!=0:
            outFile.write('C_pos = '+str(c_pos)+'\n')
            outFile.write('C_neg = '+str(c_neg)+'\n')
        #--------------------------------------------------------------
        # Plot results
        #--------------------------------------------------------------
        th=th-th0_mid
        thint=thint-th0_mid
        if opts.plotting == 1:
            fig = plt.figure(count)
            #plt.title('Peak')
            ax1=fig.add_subplot(111)
            ###########################################################
            # assign labels
            if nder > 0:
                label0='der. '+fn
            else:
                label0=fn
            ###########################################################
            ax1.plot(cf*th,r,'bo',label=label0)
            ax1.plot(cf*thint,r_g, color='r', label='Gaussian')
            ax1.plot(cf*thint,r_l, color='g', label='Lorentzian')
            ax1.plot(cf*th,bkg+th*0, color='k', label='background')
            ax1.plot(cf*th,0.5*(r_max+bkg)+th*0, linestyle='--', color='k', label='half max.')
            #
            ax1.legend(loc='upper left')
            ax1.set_xlabel(opts.nx, size=20)
            ax1.set_ylabel(opts.ny, size=20)
            ###########################################################
            # plot original data if derivative was taken
            if nder > 0:
                ax2=ax1.twinx()
                thplot=cf*(th_orig-th_max)
                rplot=r_orig
                ax2.plot(thplot,rplot,color='m',label=fn)
                #
                ax1.set_ylabel('Derivative, [a.u.]', size=20)
                ax2.legend(loc='upper right')
                ax2.set_ylabel(opts.ny,color='m', size=20)
                for tl in ax2.get_yticklabels():
                    tl.set_color('m')
                ax2.yaxis.set_major_formatter(ScalarFormatter(useOffset=True)) #(FormatStrFormatter('%.1e'))
             ##########################################################
#                errorbar(x,y,utexp,fmt='o', color=cbank[count],label=fn)
#                xdeb=arange(5,300,0.2)
#                ydeb=t3v(p,xdeb)
#                plot(xdeb,ydeb, linestyle='--', color=cbank[count], label=fn+' empirical fit function')
#                axis([0, 100, -1e-8, 0.5*1e-7])
            ###########################################################
        if opts.norm==1:
            r=r/max(r)
        figure(101)
        plot(cf*th,r,label=fn)
        figure(102)
        semilogy(cf*th,r,label=fn)
        axis([min(cf*th),max(cf*th),1e-6,1.1])
#        figure(667)
#        plot(0.024*array(range(0,count)),20.0*array(thmaxg_),'o',label='peak position')
#        plot(0.024*array(range(0,count)),20.0*array(fwhmg_),'o',label='FWHM')
#        xlabel('Energy, meV')
#        ylabel('peak position, um')

    fwhm_=cf*array(fwhm_)
    fwhm_mean=mean(fwhm_)
    fwhm_std=std(fwhm_)

    com_=cf*array(com_)
    com_mean=mean(com_)
    com_std=std(com_)

    outFile.write('###################################################################\n')
    outFile.write('fwhm mean = '+str(fwhm_mean)+'\n')
    outFile.write('fwhm std  = '+str(fwhm_std)+'\n')
    outFile.write('com mean = '+str(com_mean)+'\n')
    outFile.write('com std  = '+str(com_std)+'\n')
    outFile.close()

    figure(101)
    xlabel(opts.nx)
    if nder > 0:
        ylabel('Derivative, [raw u.]')
    else:
        ylabel(opts.ny)
    title('Common plot')
    legend(loc='upper right')


    show()


if __name__ == '__main__':
    main()