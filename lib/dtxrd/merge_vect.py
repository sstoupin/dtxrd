#####################################################################
## Purpose: tools to merge and sort lists which have identical entries
##
## Author: Stanislav Stoupin <sstoupin@aps.anl.gov>
##
## Copyright 2012 Argonne National Laboratory
##
## See the file "LICENSE" for information on usage and redistribution
##  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#######################################################################################

#####################################################################
from numpy import *

#def sort_lpair(lista,listb):
  
def sortar_pair(ar1,ar2):
  ind=ar1.argsort()
  return ar1[ind], ar2[ind]    

def mergel_same(lista,listb):
  d={}
  for a,b in zip(lista,listb):
    d.setdefault(a, []).append(b)
  
  avg=[]
  for key in d:
    avg.append(sum(d[key])/len(d[key]))
  
  return d.keys(), avg

