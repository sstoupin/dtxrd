###################################################################
# area detector instrument file
# includes basic instrument parameters 
# and signal conditions
# Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
###################################################################
xsize = 13.312      # CCD horizontal sensor size
ysize = 13.312     # CCD vertical sensor size
mfactor = 1     # optical magification factor
Mx=1024           # nominal number of x pixels
My=1024           # nominal number of y pixels
rbin=4
#
dx = xsize/mfactor/float(Mx)  # x pixel size [mm]
dy = ysize/mfactor/float(My)  # y pixel size [mm]
#
tot_range=100.0   # upper limit for threshold processing (to exclude "dead" pixels)
dyn_range=100.0   # upper limit for threshold processing (gauss)
fwhm_0=25         # [urad] expected fwhm (for fitting)
bkg0=500          # "reasonable" background - dark current count
#
data_path='/entry/instrument/detector/data'
th_path='/entry/instrument/NDAttributes/theta'
chi_path='/entry/instrument/NDAttributes/chi' 

#th0 = 428.96271936
#a2urad = 1.0
