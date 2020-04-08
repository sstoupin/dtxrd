###################################################################
# area detector instrument file
# includes basic instrument parameters 
# and signal conditions
# Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
###################################################################
xsize = 13.312      # CCD horizontal sensor size
ysize = 13.312     # CCD vertical sensor size
mfactor = 1     # optical magification factor
Mx=1024/4           # nominal number of x pixels
My=1024/4           # nominal number of y pixels

# Paths for reading data from beamline .h5 files
data_path='/entry/instrument/detector/data'      # default for 1-BM
th_path='/entry/instrument/NDAttributes/theta'   # default for 1-BM
chi_path='/entry/instrument/NDAttributes/chi'    # default for 1-BM

# Parameters used calculation of topographs
rbin=1   # NxN binning factor 
dx = xsize/mfactor/float(Mx)  # x pixel size [mm]
dy = ysize/mfactor/float(My)  # y pixel size [mm]
uf = 'urad'     # units: deg, urad, arcsec (default: deg)
zf = 1          # card for representation of intensity topographs 
#              (0 peak value normalized by max, 1 - raw detector counts, -1 integrated and normalized by max)
rng = '1.0 12.5 4.8 8.8'  # range for analysis 'x1 x2 y1 y2' in [mm]
trans = 1                 # transpose? 1 -yes 0 - no 
stat = None

# Parameters for presentation (plotting) of topographs
sample_name = 'Diamond 1'
tr = 10.0         # lower intensity rejection threshold (signal relative to background) for plotting
dyn_range = 1e3   # upper rejection threshold (rejection occurs at dyn_range*bkg0)
sf = 0.1          # scale factor for colormap range on FWHM and STDEV topographs
mf = 0.1          # scale factor for colormap range on COM, Midpoint, Left Slope and Right Slope topographs
fwhm_0=25         # [urad] expected fwhm (for fitting) (redundant, not required by rctopo-fast)
bkg0=100          # "reasonable" background - dark current count
