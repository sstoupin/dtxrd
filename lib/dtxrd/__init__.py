# Name: dtxrd/__init__.py
# Purpose: Main module for the DTXRD package
# Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
#
# Copyright 2014
#
# See the file "LICENSE" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.

"""
The main module for DTXRD.

@version: 1.0.0

@author: Stanislav Stoupin

@contact: sstoupin@aps.anl.gov

@copyright: 2014

@license: MIT X11/XFree86 style (see the file C{LICENSE} for more information)
"""

__version__ = '1.0.0'

from numpy import *
from scipy import *
from scipy.special import legendre
from pylab import *

from myio import *
from curvestat import *
from thfind import *
from dtxrd import *
from constants import *
from chi import *
from rotation import *
from dtxrd2_k import *

from pyhdf4 import *
from fit1d import *
from matplotlib.patches import Patch

#import dtxrd.pyhdf4 as pyhdf4
#import dtxrd.myio as myio
