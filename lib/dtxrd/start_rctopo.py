#!/usr/bin/env python

'''
a program to calculate rocking curve images for a series of hdf4 snapshots

:author:    Stanislav Stoupin
:email:     sstoupin@aps.anl.gov

:copyright: Copyright 2014 by XSD, Advanced Photon Source, Argonne National Laboratory
:license:   UChicago Argonne, LLC Open Source License, see LICENSE for details.
'''

import sys
import numpy
import operator
from pylab import *

from pyhdf4 import *
from curvestat import *
from fit1d import *
from myio import *

#from __future__ import division
from matplotlib.patches import Patch
##################################################################
## constants
##################################################################
r2d=180.0/pi

dx, dy = 0.06, 0.06  # CCD camera pixel size [mm]
tot_range=10.0         # upper limit for threshold processing (to exclude "dead" pixels)
dyn_range=6.0      # upper limit for threshold processing (gauss)
fwhm_0=30.0          # expected fwhm (for fitting)
fwhm_l=0.5          # lower fract. limit to plot fwhm distrib
fwhm_r=1.5           # upper limit to plot fwhm distribution
bkg0=2750.0          # "reasonable" background - dark current count

##################################################################

__version__='0.67'
proginfo = 'topohdf' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>'

# version 0.2 - added gaussian fit for noisy data as an option "-g"

#######################################################################################
# HISTORY
#######################################################################################
# 1. need deglitching algorithm - done v 0.4
# 2. -b - custom background done v 0.4
# 3. calculate RMS and PV - done in v0.3
# 4. convenient output  done in v 0.3
# 5. a bug fixed with image orientation v 0.4 - for some reason the arrays have to be
# transposed  v0.4
# 6. introduced colormap range factor (-f _value_  option)
# 7. statistics for selected rectangular region implemented   done v 0.6
# 9. need to sort data for every pixel based on angle values! done v 0.5
# 10. added upper limit on dyn_range (fixed value)            done v 0.51
# 11. added search for peak index to improve gaussian fit     done v 0.52
# 12. move dyn_range and other constants to the beginning
#       introduce df as an input parameter to the parser      done v 0.61
#       deglitch by default                                   done v 0.63
#       improve threshold upper limit processing              done v 0.63
#       subtract bkg from peak data                           done v 0.63
#       add bkg0 - dark current                               done v 0.64
# 13. improve glitch and background determination (4 points)  done v 0.65
# v 066  fliplr applied to flip left to right so that the images show up as the
#        incident beam sees them
#        added sample name to report plot as an option
# v067   additional pixel rejection criteria based on the bkg level and on the
#        COM and FWHM values (worked with Si 220 transmission data)
#        eventually best way to reject bad pixels - just treshold (e.g. -t 1.2)
#######################################################################################
#  TODO LIST
#######################################################################################
# 0. need to improve comput. speed!
# 1. need mouse click access to each point on the maps
# 2. need total rocking curve
# 3. need to select region for processing
#######################################################################################

def residuals(p,y,x):
    return y-gauss(p,x)

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

    USAGE = '%prog [OPTIONS...] filename1 filename2 ... \n' \
    +'the program computes and displays topography images using multiple CCD snapshots\n' \
    +'of reflectivity (hdf files) taken at different angles on the crystal rocking curve\n' \
    +'options:\n' \
    +'-o - output file; columns: x, y, peak, fwhm, com, thmid, th_neg, th_pos\n' \
    +'-t - threshold in reflectivity used to define crystal boundaries'
    VERSION = '%prog ' + __version__ + ', by Stanislav Stoupin <sstoupin@aps.anl.gov>\n'
    parser = OptionParser(usage=USAGE, version=VERSION)
    #
    # Behavior
    #

    #Output
    parser.add_option('-o', '--output',
            action='store',
            dest='output',
            default=None,
            help='write results to file F (defaults to stdout)',
            metavar='FILE')
    parser.add_option('-t', '--threshold',
            action='store',
            dest='tr',
            default=1.05,
            help='processing threshold in reflectivity to define the crystal boundaries',
            metavar='CONST')
    parser.add_option('-b', '--background',
            action='store',
            dest='bkg',
            default=-1,
            help='background from CCD',
            metavar='CONST')
    parser.add_option('-r', '--range',
            action='store',
            dest='rng',
            default='-1 -1 -1 -1',
            help='map range for analysis',
            metavar='STRING')
    parser.add_option('-n', '--name',
            action='store',
            dest='nm',
            default='',
            help='sample name',
            metavar='STRING')
    parser.add_option('-f', '--factor',
            action='store',
            dest='f',
            default=1.0,
            help='factor*fwhm_aver defines colormap range of images',
            metavar='CONST')
    parser.add_option('-d', '--deglitch',
            action='store',
            dest='degl',
            default=-1.0,
            help='deglitch data',
            metavar='CONST')
    parser.add_option('-g', '--gaussian',
            action='store_const',
            const=1,
            dest='stat',
            default=0,
            help='do gaussian curve fitting - for noisy data')
    parser.add_option('-s', '--trasnpose',
            action='store_const',
            const=1,
            dest='transpose',
            default=0,
            help='transpose images for plotting')
    parser.add_option('-u', '--units',
            action='store',
            dest='units',
            default='deg',
            help='angular units: deg, arcsec, urad',
            metavar='STRING')
    parser.add_option('-p', '--publish',
            action='store_const',
            const=1,
            dest='publish',
            default=0,
            help='generate figures for publication')

    #Parameters

    opts, args = parser.parse_args(args)
    if len(args) < 1:
        parser.print_usage()
        sys.exit(1)
    return opts, args
###############################################################################


def main():
    opts, args = ParseArguments(sys.argv[1:])
    if opts.output is not None:
        try:
            outFile = open(opts.output, 'w')
        except IOError, e:
            fatalIOError(e)
    else:
        outFile = open('topohdf.log', 'w') #sys.stdout

    outFile.write(str(opts)+' '+str(args)+'\n')
    outFile.write('topohdf version: '+__version__+'\n')
    outFile.write('range (x1 x2 y1 y2): '+opts.rng+'\n')

#-----------------------------------------------------
# initialize & get constants
#-----------------------------------------------------
    alldata=[]
    data=[]
    th_deg=[]
    tr=float(opts.tr);   print "Threshold (tr) = ", tr
    df=float(opts.degl); print "Deglitching factor (df) = ", df
    sf=float(opts.f);    print "Limit scale factor for plotting (sf) = ", sf
    nf=len(args);        print "Number of files to analyze = ", nf
#-----------------------------------------------------
    for fn in args:
        angle, size, im = read_hdf4(fn)
        nx=size[0] ;print nx
        ny=size[1] ;print ny
#  print im.dtype
        print fn+' ', angle

### sorting data by angle #######################
        data_fn=[angle,im]
        alldata=alldata+[data_fn]

    alldata=sorted(alldata,key=operator.itemgetter(0))
#################################################

    for x in alldata:
        th_deg=th_deg+[x[0]]
        data=data+[x[1]]

    data=array(data); #print data

    if opts.units=='deg':
        th=1.0e6*array(th_deg)/r2d # th_deg in degrees
    elif opts.units=='arcsec':
        th=1.0e6*array(th_deg)/(60.0)**2.0/r2d     #th_deg in arcseconds
    elif opts.units=='urad':
        th=array(th_deg)
    elif opts.units=='um':
        th=5.86*array(th_deg)

    th_beg=th[0]
    th_end=th[len(th)-1]
#  th0=0.5*(th[0]+th[len(th)-1])
    th0=0.5*(th_beg+th_end)
    th0_deg=1.0e-6*r2d*th0

####################################################
## Picking range
####################################################

    step=1

    rangexy=opts.rng.split(' ')
    rngx1=float(rangexy[0]); rngy1=float(rangexy[2])
    rngx2=float(rangexy[1]); rngy2=float(rangexy[3])

    if rngx1<0:
        indx1=0
    else:
        indx1=int(rngx1/dx)
    if rngy1<0:
        indy1=0
    else:
        indy1=int(rngy1/dy)
#---------------------------------------
    if rngx2<0 or rngx2>nx*dx:
        indx2=nx
    else:
        indx2=int(rngx2/dx)
    if rngy2<0 or rngy2>ny*dy:
        indy2=ny
    else:
        indy2=int(rngy2/dy)

#  xrange=arange(indx1,indx2,step) #;print xrange
#  yrange=arange(indy1,indy2,step)
    xrange=arange(nx-indx2,nx-indx1,step)
    yrange=arange(ny-indy2,ny-indy1,step) #;print yrange
#-----------------------------------------------------
#initialize arrays
#-----------------------------------------------------
    peak=numpy.empty((nx,ny))  #1
    fwhm=numpy.empty((nx,ny))  #2
    thmid=numpy.empty((nx,ny)) #3
    com=numpy.empty((nx,ny))   #4
    thneg=numpy.empty((nx,ny)) #5
    thpos=numpy.empty((nx,ny)) #6
    fwhm1=[]; peak1=[]; thmid1=[]; com1=[]; thpos1=[]; thneg1=[] # arrays to calclulate rms and pv
#-------------------------------------------------------
    if opts.bkg>0:
        bkg=float(opts.bkg); print "bkg = ", bkg
    else:
        print "bkg approximated by the endpoints"
#-------------------------------------------------------
    abs_max=data.max(); print 'abs_max = ', abs_max
#-------------------------------------------------------
# Main stuff
#-------------------------------------------------------
    for j in yrange:
        for i in xrange:
            rcurve=data[:,j,i]
            #-------------------------------------------------
            # Background determination
            #-------------------------------------------------
            if opts.bkg==-1:
                bkg=0.25*(rcurve[0]+rcurve[1]+rcurve[len(rcurve)-2]+rcurve[len(rcurve)-1]) #; print 'bkg = ', bkg

            #------------------------------------------------
            # Deglitching procedure
            #------------------------------------------------
            if df>1.0:
                index=range(len(rcurve))
                index=index[2:len(rcurve)-2]
                for k in index:
                    y=rcurve[k]
                    mean=0.25*(rcurve[k-2]+rcurve[k-1]+rcurve[k+1]+rcurve[k+2])
                    if y > df*mean or y < (2.0-df)*mean:
                        rcurve[k]=mean
                        print "X= ", dx*i, "Y= ",dy*j, " glitch found"
            #------------------------------------------------
            rcmax=max(rcurve)
            rmax=rcmax-bkg                   #; print 'rmax= ', rmax
            #------------------------------------------------
#     tot=sum(rcurve)/len(rcurve)      #; print 'tot = ', tot
            rcurve_l=list(rcurve)
            max_i=rcurve_l.index(rcmax)
            thi0=th[max_i]-th0

#      if rmax==abs_max:
#      if (i*dx==1.92) and (j*dy==1.92):
#     specific point chosen:
#      xi=60; yi=40
#      thi0=-7.0
            indx0=int(0.5*(xrange[0]+xrange[len(xrange)-1]))
#      int(0.5*len(xrange)) #; print indx0
            indy0=int(0.5*(yrange[0]+yrange[len(yrange)-1]))
            if (i==indx0) and (j==indy0):
                print "thi0 = ", thi0
                ugol=th-th0 #th_deg
                #fwhm_0=10.0  #1e-6*50.0*r2d
                f1=plt.figure()
                plt.plot(ugol,rcurve)
                p0=[bkg,rcmax,thi0,fwhm_0] #bkg max center std
                p, unc, dte, stat, cov = fit1d(residuals,p0,ugol,rcurve); print "p @ X0,Y0 = ", p
                com0=sum((rcurve-p[0])*(ugol))/sum(rcurve-p[0]); print "com @ X0,Y0 = ", com0
                print "FWHM @ X0,Y0 =", 2.0*abs(p[3])*sqrt(log(2)), " urad"
                plt.plot(ugol,gauss(p,ugol))
                plt.title ('pixel location [mm]: X ='+str((nx-i)*dx)+', Y = '+str((ny-j)*dy))
                plt.xlabel('angle, [urad]', size='x-large')
                plt.ylabel('Counts [a.u.]', size='x-large')
                ax = plt.gca()
                fontsize=16
                for tick in ax.xaxis.get_major_ticks():
                    tick.label1.set_fontsize(fontsize)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label1.set_fontsize(fontsize)

                if opts.output is not None:
                    header="# RC from the selected pixel \n# columns: ugol rcurve \n"
                    writeFile("1pix.dat", header, ugol, rcurve)
##################################################################################################################
####### PIXEL REJECTION CRITERIA #################################################################################
##################################################################################################################
            if bkg<dyn_range*bkg0 and rcmax>tr*bkg and rcmax<tot_range*bkg0:

                if opts.stat==1:
                    p0=[bkg,rcmax,thi0,fwhm_0]
                    p, unc, dte, stat, cov = fit1d(residuals,p0,th-th0,rcurve)
                    if abs(p[0]) < dyn_range*bkg0 and abs(p[2]) < 0.5*abs(th_end-th_beg) and abs(p[3]) < abs(th_end-th_beg):  # and cov !='NA' # as an option
        #
                        fwhm_g=2.0*abs(p[3])*sqrt(log(2))
                        fwhm[i,j]=fwhm_g;      fwhm1=fwhm1+[fwhm[i,j]]
                        peak[i,j]=p[1];        peak1=peak1+[peak[i,j]]   # p[1] - in gauss function background already subtracted
                        thmid[i,j]=p[2];       thmid1=thmid1+[thmid[i,j]]
                        com[i,j]=sum((rcurve-p[0])*(th-th0))/sum(rcurve-p[0]); com1=com1+[com[i,j]]
                        thneg[i,j]=p[2]-fwhm_g/2.0; thneg1=thneg1+[thneg[i,j]]
                        thpos[i,j]=p[2]+fwhm_g/2.0; thpos1=thpos1+[thpos[i,j]]
                    else:
                        fwhm[i,j]=1e9
                        peak[i,j]=0
                        thmid[i,j]=1e9
                        com[i,j]=1e9
                        thneg[i,j]=1e9
                        thpos[i,j]=1e9
                else:
                    try:
                        stat=curvestat(th-th0,rcurve,bkg)
                        if abs(stat[1]-bkg) < dyn_range*bkg and abs(stat[6]) < 0.5*abs(th_end-th_beg) and abs(stat[5]) < abs(th_end-th_beg):
                            fwhm[i,j]=stat[5];       fwhm1=fwhm1+[fwhm[i,j]]
                            peak[i,j]=stat[1]-bkg;   peak1=peak1+[peak[i,j]]
                            thmid[i,j]=stat[4];      thmid1=thmid1+[thmid[i,j]]
                            com[i,j]=stat[6];        com1=com1+[com[i,j]]
                            thneg[i,j]=stat[2];      thneg1=thneg1+[thneg[i,j]]
                            thpos[i,j]=stat[3];      thpos1=thpos1+[thpos[i,j]]
                        else:
                            fwhm[i,j]=1e9
                            peak[i,j]=0
                            thmid[i,j]=1e9
                            com[i,j]=1e9
                            thneg[i,j]=1e9
                            thpos[i,j]=1e9
####################################################################
## Testing
####################################################################
#                if com[i,j] > 0:
#                  print "flag!", i, j
                            #com[i,j]=-1e3
#                if (i==xi) and (j==yi):
#                    print "fwhm check", fwhm[i,j]
#                    print "com check", com[i,j]
####################################################################

                    except UnboundLocalError:
                        fwhm[i,j]=1e9
                        peak[i,j]=0
                        thmid[i,j]=1e9
                        com[i,j]=1e9
                        thneg[i,j]=1e9
                        thpos[i,j]=1e9

            else:
                fwhm[i,j]=1e3
                peak[i,j]=0
                thmid[i,j]=1e3
                com[i,j]=1e3
                thneg[i,j]=1e3
                thpos[i,j]=1e3

#  norm=peak.max()
    norm=max(peak1)
    peak=peak/norm
    peak1=peak1/norm

    fwhm1=array(fwhm1)
    peak1=array(peak1)
    thmid1=array(thmid1)
    com1=array(com1)
    thneg1=array(thneg1)
    thpos1=array(thpos1)

    print '-----------------------------------------------------------------------'
    out=''
    #----------------------------------------------------------------------------------
    fwhm0=average(fwhm1);          out=out+"FWHM_average = "+str(fwhm0)+" urad\n"
    fwhm_rms=std(fwhm1);           out=out+"FWHM_rms = "+str(fwhm_rms)+" urad\n"
    fwhm_pv=max(fwhm1)-min(fwhm1); out=out+"FWHM_pv = "+str(fwhm_pv)+" urad\n"
    #----------------------------------------------------------------------------------
    peak0=average(peak1);          out=out+"Reflect_average = "+str(peak0)+"\n"
    peak_rms=std(peak1);           out=out+"Reflect_rms = "+str(peak_rms)+"\n"
    peak_pv=max(peak1)-min(peak1); out=out+"reflect_pv = "+str(peak_pv)+"\n"
    #----------------------------------------------------------------------------------
    com0=average(com1);         out=out+"COM_average = "+str(com0)+" urad\n"
    com_rms=std(com1);          out=out+"COM_rms = "+str(com_rms)+" urad\n"
    com_pv=max(com1)-min(com1); out=out+"COM_pv = "+str(com_pv)+" urad\n"
    #----------------------------------------------------------------------------------
    thmid0=average(thmid1);           out=out+"thmid_average ="+str(thmid0)+" urad\n"
    thmid_rms=std(thmid1);            out=out+"thmid_rms = "+str(thmid_rms)+" urad\n"
    thmid_pv=max(thmid1)-min(thmid1); out=out+"thmid_pv = "+str(thmid_pv)+" urad\n"
    #----------------------------------------------------------------------------------
    thneg0=average(thneg1);           out=out+"thneg_average = "+str(thneg0)+" urad\n"
    thneg_rms=std(thneg1);            out=out+"thneg_rms = "+str(thneg_rms)+" urad\n"
    thneg_pv=max(thneg1)-min(thneg1); out=out+"thneg_pv = "+str(thneg_pv)+" urad\n"
    #----------------------------------------------------------------------------------
    thpos0=average(thpos1);           out=out+"thpos_average = "+str(thpos0)+" urad\n"
    thpos_rms=std(thpos1);            out=out+"thpos_rms = "+str(thpos_rms)+" urad\n"
    thpos_pv=max(thpos1)-min(thpos1); out=out+"thpos_pv = "+str(thpos_pv)+" urad\n"
    #----------------------------------------------------------------------------------
    print out
##################################################################################
## Write output                                                                 ##
##################################################################################

    outFile.write('-----------------Statistics---------------------------------------------------------\n')
    outFile.write(out)
    outFile.write('-----------------Map of the selected range------------------------------------------\n')

    for j in yrange:
        for i in xrange:
            outFile.write(str(i*dx)+' '+str(j*dy)+' '+str(peak[i,j])+' '+str(fwhm[i,j])+' ' \
            +str(com[i,j])+' '+str(thmid[i,j])+' '+str(thneg[i,j])+' '+str(thpos[i,j])+'\n')

    outFile.close
    print 'peak.max = ', norm
    print tr
##################################################################################
##  Plot results
##################################################################################
    # need to "mirror" y indexes because in an array numbering starts from top left corner
    # affects only choice of range - does not affect plotting
    # same for x since it was found that the image has to be reversed in x direction
    # note that for x fliplr is applied later while for y imshow origin='upper'

    nny1=ny-indy2; nny2=ny-indy1
#  nny1=indy1; nny2=indy2
    nnx1=nx-indx2; nnx2=nx-indx1
#  nnx1=indx1; nnx2=indx2

    fwhm=fwhm[nnx1:nnx2,nny1:nny2]
    peak=peak[nnx1:nnx2,nny1:nny2]
    thmid=thmid[nnx1:nnx2,nny1:nny2]
    com=com[nnx1:nnx2,nny1:nny2]
    thneg=thneg[nnx1:nnx2,nny1:nny2]
    thpos=thpos[nnx1:nnx2,nny1:nny2]


    if opts.transpose==1:
        fwhm=transpose(fwhm)
        peak=transpose(peak)
        thmid=transpose(thmid)
        com=transpose(com)
        thneg=transpose(thneg)
        thpos=transpose(thpos)
#    xyrange=(0,dx*nx,0,dy*ny)  # for some reason the arrays should be transposed
#  else:
#    xyrange=(0,dy*ny,0,dx*nx)
        xyrange=(dx*indx1,dx*indx2,dy*indy1,dy*indy2)  # for some reason the arrays should be transposed
    else:
        xyrange=(dy*indy1,dy*indy2,dx*indx1,dx*indx2)

    fwhm=fliplr(fwhm)
    peak=fliplr(peak)
    thmid=fliplr(thmid)
    com=fliplr(com)
    thneg=fliplr(thneg)
    thpos=fliplr(thpos)


    f2=plt.figure()
    if opts.stat==1:
        sbtitle='Gaussian curve fit.'
    else:
        sbtitle='No fitting.'

    sbtitle='Rocking curve topographs: '+sbtitle

    if opts.nm != '':
        sbtitle=sbtitle+' Sample: '+opts.nm

    plt.suptitle(sbtitle, size=20)

    plt.subplot(331)
    plt.imshow(peak, aspect='auto', extent=xyrange,vmin=0,vmax=1.0)
    plt.colorbar()
    plt.title('normalized reflectivity',size=20)
    #plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

#  vmin0=round(fwhm_l*fwhm0,-1)
#  vmax0=round(fwhm_r*fwhm0,-1)
    vmin0=round(fwhm0*(1.0-sf))
    vmax0=round(fwhm0*(1.0+sf))

    plt.subplot(332)
    plt.imshow(fwhm, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=fwhm.min(),vmax=fwhm.max())
    plt.colorbar()
    plt.title('FWHM [urad]',size=20)
    #plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    #grid(True)

    vmin0=round(com0-sf*fwhm0)
    vmax0=round(com0+sf*fwhm0)

    plt.subplot(334)
    plt.imshow(com, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=com.min(),vmax=com.max())
    plt.colorbar()
    plt.title('Center of mass (angle) [urad]',size=20)
    #plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    vmin0=round(thmid0-sf*fwhm0)
    vmax0=round(thmid0+sf*fwhm0)

    plt.subplot(335)
    plt.imshow(thmid, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=thmid.min(),vmax=thmid.max())
    plt.colorbar()
    plt.title('mid-point (angle) [urad]',size=20)
    #plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    vmin0=round(thneg0-1.0*sf*fwhm0)
    vmax0=round(thmid0)

    plt.subplot(337)
    plt.imshow(thneg, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=thneg.min(),vmax=thneg.max())
    plt.colorbar()
    plt.title('left slope (angle) [urad]',size=20)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

#  vmin0=round(com0-(sf-1.0)*fwhm0,-1)
    vmin0=round(thmid0)
    vmax0=round(thpos0+1.0*sf*fwhm0)

    plt.subplot(338)
    plt.imshow(thpos, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=thpos.min(),vmax=thpos.max())
    plt.colorbar()
    plt.title('right slope (angle) [urad]',size=20)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')


    plt.subplot(333)
    axis('off')
    plt.title('X-ray beam orientation - diffraction in horizontal plane')
    plt.annotate('vertical', xy=(0.5,1.0), xytext=(0.5,0.0),
                arrowprops=dict(facecolor='black', shrink=0.05),
                            )
    plt.annotate('horizontal', xy=(1.0,0.5), xytext=(0.0,0.5),
                arrowprops=dict(facecolor='black', shrink=0.05),
                            )
    plt.subplot(339)
    axis('off')
    plt.text(0,0,out,size=20)
#  plt.show()
    ################################
    # output for GLE
    ################################
#  thc=59.84
#  dx1=dx/sin(thc/r2d)
    if opts.output is not None:
        com_exp=array(com) #-com0
        #
        outfile1 = open('com.dat', 'w')
        outfile2 = open('fwhm.dat', 'w')
        outfile3 = open('com1.dat', 'w')
        scom=shape(com_exp) #print "scom =", scom
#    sfwhm=shape(fwhm)
        for i in range(0,scom[0],1):
            for j in range(0,scom[1],1):
                outfile1.write(str(i*dx)+' '+str(j*dy)+' '+str(com_exp[i,scom[1]-1-j])+'\n')
                outfile2.write(str(i*dx)+' '+str(j*dy)+' '+str(fwhm[i,scom[1]-1-j])+'\n')

        for x in com1:
            outfile3.write(str(x)+'\n')
        outfile1.close
        outfile2.close
        outfile3.close
    if opts.publish==1:
        from publish import figplot
        com_exp=array(com)-com0
        figplot(opts.transpose,dx,dy,indx1,indx2,indy1,indy2,com_exp,fwhm)
    plt.show()
#
#    matplotlib.rcParams.update({'font.size': 20})
#    f3=plt.figure()
#    vmin0=round(com0-sf*fwhm0,-1); print vmin0
#    vmax0=round(com0+sf*fwhm0,-1); print vmax0
#    vmin0=-25
#    vmax0=25
        #
#    imgplot = plt.imshow(com, aspect='auto', extent=xyrange, vmin=vmin0, vmax=vmax0) #vmin=com.min(),vmax=com.max())
#    imgplot.set_cmap('spectral')
#    plt.colorbar()
        #plt.title('Center of mass (angle) [urad]',size=20)
#    plt.xlabel('x [mm]')
#    plt.yticks([1,2,3,4,5],('1.0','2.0','3.0','4.0','5.0'))
#    plt.ylabel('y [mm]')
#    plt.ylim(ymin=0)
#    plt.spines['left'].set_position('zero')
#
#  plt.show()


if __name__ == '__main__':
    main()