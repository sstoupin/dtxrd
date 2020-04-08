###################################################################
# area detector instrument file
# includes basic instrument parameters 
# and signal conditions
# Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
###################################################################
xsize = 16.6      # CCD horizontal sensor size
ysize = 14.0      # CCD vertical sensor size
mfactor = 10.0     # optical magification factor
Mx=2560/4           # nominal number of x pixels
My=2160/4           # nominal number of y pixels
rbin=1
#
#dx = xsize/mfactor/float(Mx)  # x pixel size [mm]
#dy = ysize/mfactor/float(My)  # y pixel size [mm]
dx = 0.00259375
dy = 0.00259375
#
tot_range=1000.0       # upper limit for threshold processing (to exclude "dead" pixels)
dyn_range=1000.0        # upper limit for threshold processing (gauss)
fwhm_0=3.0          # [urad] expected fwhm (for fitting)
bkg0=120.0          # "reasonable" background - dark current count
#
data_path='/entry/instrument/detector/data'
th_path='/entry/instrument/NDAttributes/theta'
chi_path='/entry/instrument/NDAttributes/chi' 
