#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:26:00 2025

BacDraw reads binary data written by FrameAnalysisTool (=FAT) and generates a plot.
This may be used as a template for an interface with FAT data.

@author: jhumt


 This file is part of FrameAnalysisTool
 

     FrameAnalysisTool, is free  software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the 
     Free Software Foundation, either version 3 of the License or any later 
     version.
 
     FrameAnalysisTool is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License along 
     with FrameAnalysisTool.  If not, see <http://www.gnu.org/licenses/>.

 Copyright 2025 Juergen Humt
 
"""
import _pickle as pickle
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from utilPY import biarc_plot


def clcDst(vX, vR, pX):
# distance between a point and a straight line
# vX  point through which line passes
# vR  derection of line
# pX  point in question
#
  vP = pX - vX
  d = abs(vP.dot(vR))
  vY = vP - d*vR
  dT= norm(vY)
  return dT


wkDir = '/home/jhumt/Flight/FreeCad/Macro_1/FrameAnalysisTool/'

    
inpFname = 'frm6_strm'


bacFname = wkDir + inpFname + '_bac.dat'

firstrunDICT = True
if firstrunDICT:
  bacF = open(bacFname,'rb')
  bacData = pickle.load(bacF)
  bacF.close()



plt.figure(figsize=(8, 8))       
# plt.axhspan(-40, 170)
# plt.axvspan(-80, 10)

fig, axs = plt.subplots(1,1)
axs.set_aspect("equal")


fsg1={"vS": np.array([0,0]), "vE": np.array([0,0])}
for j in range(len(bacData) ):
  bacItm = bacData[j]
  fsg1['vS']=bacItm["fsg1vS"]
  fsg1['vE']=bacItm["fsg1vE"]
  l0 =bacItm["l0"]
  theta0 =bacItm["theta0"]
  kappa0 =bacItm["kappa0"]
  
  
  l1 =bacItm["l1"]
  theta1 =bacItm["theta1"]
  kappa1 =bacItm["kappa1"]  
  
  biarc_plot(fsg1['vS'][0], fsg1['vS'][1], l0, theta0, kappa0, fsg1['vE'][0], fsg1['vE'][1],  l1, theta1, kappa1, fmt1=None, fmt2=None)  
  
  # : , "l0": l0, "theta0": theta0, "kappa0": kappa0, "fsg1vE": fsg1['vE'], "l1": l1, "theta1": theta1, "kappa1": kappa1, "xs": xs, "ys": ys,  "thetas": thetas }
  k=0
axs.grid(color='b', linestyle='-') # b=True, which='major', color='b', linestyle='-')
plt.show()