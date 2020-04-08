###################################################################
# area detector instrument file
# includes basic instrument parameters 
# and signal conditions
# Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
###################################################################
xsize = 16.6      # CCD horizontal sensor size
ysize = 14.0      # CCD vertical sensor size
mfactor = 1.0     # optical magification factor
Mx=2560           # nominal number of x pixels
My=2160           # nominal number of y pixels


# Parameters used calculation of topographs
rbin=1   # NxN binning factor 
dx = xsize/mfactor/float(Mx)  # x pixel size [mm]
dy = ysize/mfactor/float(My)  # y pixel size [mm]
uf = 'deg'     # units: deg, urad, arcsec (default: deg)
#zf = 1          # card for representation of intensity topographs 
#              (0 peak value normalized by max, 1 - raw detector counts, -1 integrated and normalized by max)
rng = '2.5 13.1 8.6 12.0' # range for analysis 'x1 x2 y1 y2' in [mm]
trans = 1                 # transpose? 1 -yes 0 - no
fn_ang = '5.dat'

# Parameters for presentation (plotting) of topographs
sample_name = 'SiC4H-6_008'
tr = 1.3          # lower intensity rejection threshold (signal relative to background) for plotting
fwhm_0=50         # [urad] expected fwhm (for fitting) (redundant, not required by rctopo-fast)
bkg0=100          # "reasonable" background - dark current count
dyn_range = 1e3   # upper rejection threshold (rejection occurs at dyn_range*bkg0)
sf = 0.5          # scale factor for colormap range on FWHM and STDEV topographs
mf = 5.0          # scale factor for colormap range on COM, Midpoint, Left Slope and Right Slope topographs
tot_range = dyn_range

#data_path = ''
#th_path = ''
#chi_path = ''