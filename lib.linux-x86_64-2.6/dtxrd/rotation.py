#!/usr/bin/python

from numpy import *

#rotation of xyz around x-axis
def rot_x(th):
      trx_1=[1.0, 0.0,      0.0    ]
      trx_2=[0.0, cos(th), -sin(th)]
      trx_3=[0.0, sin(th),  cos(th)]
      return array([trx_1,trx_2,trx_3])

# rotation of xyz around y-axis
def rot_y(th):
      try_1=[cos(th),  0.0, sin(th)]
      try_2=[0.0,      1.0, 0.0    ]
      try_3=[-sin(th), 0.0, cos(th)]            
      return array([try_1,try_2,try_3])   

#rotation of xyz around z-axis
def rot_z(th):      
      trz_1=[cos(th), -sin(th), 0.0]
      trz_2=[sin(th),  cos(th), 0.0]
      trz_3=[ 0.0,     0.0,     1.0]
      return array([trz_1,trz_2,trz_3])       

#inversion of x1
inv_x=array([[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
inv_y=array([[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0]])
inv_z=array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0]])
identity3=array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])


      