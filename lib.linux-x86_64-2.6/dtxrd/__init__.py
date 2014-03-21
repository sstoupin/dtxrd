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

from dtxrd import myio
from myio import *
from dtxrd.curvestat import *
from dtxrd.thfind import *
from dtxrd.dtxrd import *
from dtxrd.constants import *
from dtxrd.chi import *
from dtxrd.rotation import *
from dtxrd.dtxrd2_k import *

from dtxrd.pyhdf4 import *
from dtxrd.fit1d import *
from matplotlib.patches import Patch
