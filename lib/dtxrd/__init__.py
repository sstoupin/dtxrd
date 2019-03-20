# Name: dtxrd/__init__.py
# Purpose: Main module for the DTXRD package
# Author: Stanislav Stoupin <sstoupin@gmail.com>
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

@copyright: Copyright (c) 2014 by XSD, Advanced Photon Source, Argonne National Laboratory

@license: UChicago Argonne, LLC OPEN SOURCE LICENSE
"""

__version__ = '1.0.0'

from numpy import *
from scipy import *
#from scipy.special import legendre
from pylab import *

from myio import *
from curvestat import *
from thfind import *
from dtxrd0 import *
from constants import *
from chi import *
from rotation import *
from dtxrd2_k import *
from pyhdf import *
from fit1d import *
from matplotlib.patches import Patch
