#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:59:08 2025

FrameAnalysisTool allows to calculate biarc approximations
to ship frames.

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

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, acos, atan, atan2, degrees, radians, pi
from utilPY import biarc_plot, atan36, circlefit
from BiarcM import biarc
from svg.path import parse_path
from xml.dom import minidom
from sys import exit
import _pickle as pickle





def clcSegLen(xAr, zAr):
  alfBnd = radians(5)
  epsB= 1.0 - 1.0e-5
  vS = np.array([xAr[1],zAr[1]]) - np.array([xAr[0],zAr[0]]) 
  vS = vS/np.linalg.norm(vS)    
  k=0
  k_p=[]
  # for x,z in zip(xAr,zAr):
  for j in range(1,len(xAr)-1):
     x=xAr[j+1] - xAr[j]
     z=xAr[j+1] - zAr[j]
     
     vK = np.array([x,z])
     vK = vK/np.linalg.norm(vK)
     alfK = vS.dot(vK)
     if alfK > epsB:
        alfK = 1.0
     else:   
        alfK = acos(alfK)
     if alfK > alfBnd:
       k_p.append(k)
       k=0
       vS = vK
     else:
       k+=1  
  return k_p

def getSvgItm(doc,tag,attr):
  kRet = [path.getAttribute(attr) for path
             in doc.getElementsByTagName(tag)]  
  return kRet
    
def chkSvgItm(xItm):
  if len(xItm) > 0:
       try:
          cX = float(xItm[0])
       except ValueError:
           cX = 1.0
  else:
       cX = 1.0
  return cX


def mkKeel(xP,zP,pAr,zW,plt):
# the two values in pAr are zK, the depth of the keel
# below the last frame value and xK, the width of the
# frame. If xK=0 the lower keel line will end at x=0
#

   if len(pAr[0]) > 2:
      pAr=pAr[0].split(',')   
      xK=float(pAr[0])
      zK=float(pAr[1])
      jK=int(pAr[2])

      if jK >= 0:
        jK=0

      x1 = xP[jK]
      z1 = zP[jK]
      
      x2 = x1
      z2 = zP[jK] + zK

      if abs(xK) > 1.0e-6:
        x3 = x1 + xK
      else:
        x3 = 0
        
      z3 = z2
        
      p1 = ([x1,x2],[z1,z2])
      p2 = ([x2,x3],[z2,z3])

      plt.plot([p1[0],p2[0]],[p1[1],p2[1]],'-b', linewidth=2)      



wkDir = '/home/jhumt/Flight/FreeCad/Macro_1/FrameAnalysisTool/data/'

fig, axs = plt.subplots(1,1)
size = fig.get_size_inches()*fig.dpi


# Read data, change value below (remove minus) 
# to read from csv file. 
useCsv = -1
if useCsv > 0:
    
    # 
    inpFname='frm6_strm2'
    
    # the two segmentations below give odd results
    k_p = [8, 9, 7, 6, 6]
    k_p = [8, 9, 3, 4, 6, 6]
    
    # the segementation below gives a proper result
    k_p = [4, 6, 14, 4, 5, 4]    

    
    f_name = wkDir + inpFname + '.csv'
    data = np.loadtxt(f_name, delimiter=',')
    
   
    outName = wkDir + inpFname + '_out.lis'
    outF = open(outName,'w')
    
else:   

    # inpFname = 'frm6_strm2' loop    
    inpFname='frm6_strm'
    # inpFname = 'fr4_dutch'
    inpFname = 'fr4_dutch2'
    # inpFname = 'fr6_dutch4'
    # inpFname = 'fr6_dutch4c'
    
    print('data from file: ' + inpFname)
    f_name = wkDir + inpFname + '.svg' 
    
    # read the SVG file
    doc = minidom.parse(f_name)
    path_strings = [path.getAttribute('d') for path
                    in doc.getElementsByTagName('path')]
    
    cX = getSvgItm(doc,'usr', 'cX')
    cX = chkSvgItm(cX)   
       
    cZ = getSvgItm(doc,'usr', 'cZ')
    cZ = chkSvgItm(cZ)


    xzK = getSvgItm(doc, 'usr', 'keel')
      
    # input from svg file should be something like:
    #   <usr kp="9 6 8 16 20 21 17" />     
    kp = getSvgItm(doc,'usr','kp')
    

    if len(kp) > 0:
        kpL = kp[0].split(' ')
        k_p = []
        for k in kpL:
           k_p.append(int(k))
    else:
    # try to provide a default, will probably 
    # not work very well.
        k_p = [19, 12, 22, 21, 17,4]
    

    doc.unlink()
    
    
    data2=[]
    eps= 1.0e-5; xB = 1.0e6; yB = 1.0e6
    for path_string in path_strings:
        path = parse_path(path_string)
        for e in path:
           x0 = e.start.real
           y0 = e.start.imag
           if (abs(xB-x0) > eps) and (abs(yB-y0) > eps):
             data2.append([x0,y0])
           xB=x0; yB=y0
           
           
    # to convert data2 to numpy array it has to
    # be split in two separate columns
    x_ = [ cX*x[0] for x in  data2]
    z_ = [ cZ*x[1] for x in  data2]
    
    # gimp uses a coordinate system, where the x values 
    # increase from right to left, make that left to right
    xMax = max(x_)
    for i in range(len(x_)):
      x_[i] = -(xMax-0.99999*x_[i])

   
    zMax=max(z_)
    for i in range(len(z_)):
      z_[i] = 0.9*zMax-z_[i]
    
    # unkomment for inverse sequence
    x_ = x_[::-1]
    z_ = z_[::-1]
    
    data = np.array([x_,z_]).T
    
  
    outName = wkDir + inpFname + '_out.lis'
    outF = open(outName,'w')
    
    
l_d = len(data)   # 
l_x = l_d
x_c = data[:l_d, 0]
z_c = data[:l_d, 1]

k_X=clcSegLen(x_c, z_c)

if useCsv > 0:
  c_m = 0.5 * 3810 / max(abs(x_c))
  x_c = c_m * x_c
  z_c = c_m * z_c

np_val = 4
print("  j   jS  jE       xCnt        yCnt             R          vTx         vTz")
outStr = "  j   jS  jE       xCnt        zCnt         vSx         vSz        vEx         vEz           R       v1Tx       v1Tz     alf1T   v2Tx     v2Tz    alf2T     zT\n" 
outF.write(outStr)

# the list below kontains the number of data points
# used to approximate each circle segment
j_s = 0  # Python uses 0-based indexing

sumKp=0
for k in k_p:
    sumKp+= k
# all remaining points form the last segment
lL = l_d-sumKp-1

outStr= f"number stations in curve {l_d:3d} number stations required by k_p {sumKp+3: 5d}\n"
print(outStr)
if lL < 3:
  exit('last segment has less than 3 elements, modify k_p')
k_p.append(lL)


plt.plot(x_c, z_c)
plt.plot(x_c, z_c, 'gx')

zW=0
mkKeel(x_c, z_c, xzK, zW, plt)

# in this section a circle is fitted to each 
# segment of the input data
lJe = len(x_c)


alfES3 = []
circData = []
varXY=float(0.0); prvAlfR = float(0.0); pi2 = 0.5*pi
frSeg={"r": varXY, "alfS": pi2, "alfE": float(90.0), "y": varXY, "z":varXY, "sTyp":varXY, "vC": np.array([varXY,varXY]), "vR": np.array([varXY,varXY]), "vS": np.array([varXY,varXY]), "vE": np.array([varXY,varXY]), "vA": np.array([varXY,varXY]), "vZ": np.array([varXY,varXY]), "alfR": prvAlfR, "alfM": prvAlfR, "bta2": float(0.0), "bta3": float(0.0), "dIr": 0}
# circData.append(frSeg)
for j in range(len(k_p)):
    j_e = j_s + k_p[j]
    plt.plot(x_c[j_e], z_c[j_e], 'rx')

    x = x_c[j_s:j_e+1]
    y = z_c[j_s:j_e+1]
    lx = len(x) - 1
    
    # calculate an approximate vector of the direction
    # at the start of each segment
    frSeg['vA']=np.array([x[1]-x[0],y[1]-y[0]])
    frSeg['vA']=frSeg['vA']/np.linalg.norm(frSeg['vA'])

    frSeg['vZ']=np.array([x[lx]-x[lx-1],y[lx]-y[lx-1]])
    frSeg['vZ']=frSeg['vZ']/np.linalg.norm(frSeg['vZ'])
    
    x_center, z_center, radius = circlefit(x, y)
    
    
    plt.plot(x_center, z_center, 'bx')
    plt.plot(x_center, z_center, 'ro')
    plt.text(x_center+6,z_center,str(j+1))

    v1 = [x[0] - x_center,y[0] - z_center]
    v1 = v1/np.linalg.norm(v1)
    # generate vecotor perpendicular to v1
    V1T= np.array([-v1[1], v1[0]])
    aA = acos(V1T.dot(frSeg['vA']))
    if aA > 0.5001*pi:
        V1T = -V1T
        aA = acos(V1T.dot(frSeg['vA']))
    
    alf1T = degrees(atan36(V1T[1],V1T[0]))   

    v2 = [x[lx] - x_center,y[lx] - z_center]    
    v2 = v2/np.linalg.norm(v2)
    
    V2T= np.array([-v2[1], v2[0]])

    aB = acos(V2T.dot(frSeg['vZ']))
    if aB > 0.5001*pi:
      V2T = -V2T
      aB = acos(V2T.dot(frSeg['vZ']))
      
    c_alf = np.dot(v1, v2) / (np.linalg.norm(v1))
    alfR = acos(c_alf)
    lC = radius*alfR
    alf = degrees(alfR)
    # vS, vE are the vectors from the center of the circle to the start and end 
    # datapoints, i.e they are perpendicular to the circle perimeter
    vS = np.array([x[0],y[0]]);       vE = np.array([x[lx],y[lx]])
    zT = np.cross(V1T,v2)
    # print(f"{j+1:3d}  {j_s+1:3d}  {j_e+1:3d}  {x_center:9.2f}   {z_center:9.2f}    {radius:9.2f}   {x_c[j_e]:9.2f}    {z_c[j_e]:9.2f}  {alf:6.1f}")
    print(f"{j+1:3d}  {j_s:3d}  {j_e:3d}  {x_center:9.2f}   {z_center:9.2f}    {radius:9.2f}   {x_c[j_e]:9.2f}    {z_c[j_e]:9.2f}  {alf:6.1f}  {v2[0]:10.5f}  {v2[1]:10.5f}")
    # alf2T = degrees(atan(v2[1]/v2[0]))
    alf2T = degrees(atan36(V2T[1],V2T[0]))

    alfES3.append([alf1T, alf2T])  
    
    # if jTest > 0:
    je1=j_e+1
    if je1==lJe:
      je1= j_e
    vE_ = np.array([x_c[je1], z_c[je1]])        
    frSeg={"r": radius, "alfS":alf1T, "alfE":alf2T, "y": x_center, "z":z_center, "sTyp":'tc', "vC": np.array([x_center, z_center]), "vR": np.array([varXY, varXY]), "vS": np.array([x_c[j_s], z_c[j_s]]), "vE": np.array([x_c[je1], z_c[je1]]), "vA": np.array([varXY,varXY]), "vZ": np.array([varXY,varXY]), "alfR": prvAlfR, "bta2": float(0.0), "bta3": float(0.0), "dIr": zT, "lC": lC}
    frSeg['vC']=np.array([x_center,z_center])
    circData.append(frSeg)
      
    outStr= f"{j+1:3d}, {j_s:3d}, {j_e+1:3d}, {x_center:9.2f},  {z_center:9.2f},   {x[0]:9.2f},  {y[0]:9.2f},  {x[lx]:9.2f},  {y[lx]:9.2f},  {radius:9.2f},  {V1T[0]:7.4f},   {V1T[1]:7.4f}, {alf1T:6.1f}, {V2T[0]:7.4f}, {V2T[1]:7.4f}, {alf2T:6.1f}  {zT:6.2}\n"
    outF.write(outStr)
    j_s = j_e


mE = len(k_p)
alfM3=[]
alfM3.append(alfES3[0][0])
j=0
outStr= f"{j+1:3d}, {alfM3[j]:8.2f} \n"
outF.write(outStr)
for j in range(1,mE):
   alfM3.append(0.5*(alfES3[j-1][1] + alfES3[j][0]))
   outStr= f"{j+1:3d}, {alfM3[j]:8.2f}\n"
   outF.write(outStr)   
alfM3.append(alfES3[mE-1][1])
outStr= f"{j+2:3d}, {alfM3[j+1]:8.2f}\n"
outF.write(outStr)
   
for j in range(mE-1):
   circData[j]['vE']= circData[j+1]['vS']   

kS=0; kE = mE 
jS = kS

# all data necessary to genrate a frame are collected
# and later on stored in binary form
bacItm = {"fsg1vS": np.array([varXY,varXY]), "l0": varXY, "theta0": varXY, "kappa0": varXY, "fsg1vE": np.array([varXY,varXY]), "l1": varXY, "theta1": varXY, "kappa1": varXY, "xs": varXY, "ys": varXY,  "thetas": varXY,   }
bacData=[]


outF.close()
outName = wkDir + inpFname + '_bac.lis'
oF=open(outName,'w')
outStr = "fsg1['vS'][0],fsg1['vS'][1],l0,       theta0, 1/kappa0,  fsg1['vE'][0],fsg1['vE'][1], l1,         theta1,        1/kappa1\n"
oF.write(outStr)

for j in range(kS,kE):
   fsg1=circData[jS]
  
   vP = 0.5*(fsg1['vS']+fsg1['vE'])
   plt.text(vP[0]+6,vP[1]-4,str(jS+1))
   
   [l0, theta0, kappa0, l1, theta1, kappa1, xs, ys, thetas] =  biarc(fsg1['vS'][0], fsg1['vS'][1], radians(alfM3[jS]), fsg1['vE'][0], fsg1['vE'][1], radians(alfM3[jS+1]))       
   # if not (oF == None):
   #     outStr = f"{jS:3d},{l0:9.4f}, {theta0:9.4f}, {kappa0:9.4f}, {l1:9.4f}, {theta1:9.4f}, {kappa1:9.4f}, {xs:9.4f}, {ys:9.4f}, {thetas:9.4f} '=  biarc(' , {fsg1['vS'][0]:9.4f},    {fsg1['vS'][1]:9.4f},       {alfM2[jS-1]:9.4f},               {fsg1['alfR']:9.4f}, {dAlf1[jS]:7.2f},     {fsg1['vE'][0]:9.4f},    {fsg1['vE'][1]:9.4f},  {alfM2[jS]:9.4f},     {fsg2['alfR']:9.4f},       {dAlf2[jS]:7.2f},       {fsg2['alfS']:9.4f},       {fsg2['alfE']:9.4f}"
   #     oF.write(outStr + '\n')
   biarc_plot(fsg1['vS'][0], fsg1['vS'][1], l0, theta0, kappa0, fsg1['vE'][0], fsg1['vE'][1],  l1, theta1, kappa1, fmt1=None, fmt2=None)
   outStr= f"{fsg1['vS'][0]:9.4f},  { fsg1['vS'][1]:9.4f},  { l0:9.4f},  { theta0:9.4f},  { 1/kappa0:8.2f},  { fsg1['vE'][0]:9.4f},  { fsg1['vE'][1]:9.4f},{  l1:9.4f}, { theta1:9.4f},  { 1/kappa1:8.2f}"
   oF.write(outStr + '\n')
   bacItm = {"fsg1vS": fsg1['vS'], "l0": l0, "theta0": theta0, "kappa0": kappa0, "fsg1vE": fsg1['vE'], "l1": l1, "theta1": theta1, "kappa1": kappa1, "xs": xs, "ys": ys,  "thetas": thetas }
   bacData.append(bacItm)
   jS+=1


if not None==oF:
  oF.close()
    
bacFname = wkDir + inpFname + '_bac.dat'
bacF = open(bacFname,'wb')

pickle.dump(bacData, bacF)
bacF.close()



j_e = 18
plt.plot(x_c[j_e], z_c[j_e], 'go')
plt.plot(x_c[j_e], z_c[j_e], 'yx')
print(f"{j_e:3d}   {x_c[j_e]:9.2f}    {z_c[j_e]:9.2f}")

j_e = 0
print(f"{j_e:3d}   {x_c[j_e]:9.2f}    {z_c[j_e]:9.2f}")

j_e = len(x_c) - 1
print(f"{j_e:3d}   {x_c[j_e]:9.2f}    {z_c[j_e]:9.2f}")



xLim = axs.get_xlim()
yLim = axs.get_ylim()


plt.text(0.7*xLim[1], 0.96*yLim[0], inpFname)
axs.set_aspect("equal")
axs.grid(color='b', linestyle='-') 
plt.show()
k = 0
