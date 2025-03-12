#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 07:22:04 2025

@author: jhumt

#  The solution algorithm is described in
#  A Note on Robust Biarc Computation
#  by Enrico Bertolazzi and Marco Frego
#  Computer-Aided Design & Applications, 16(5), 2019, 822-835
#  https://doi.org/10.14733/cadaps.2019.822-835
#  https://cad-journal.net/files/vol_16/CAD_16(5)_2019_822-835.pdf

"""

import numpy as np
from utilPY import biarc_plot
from math import degrees, sin, cos, pi
import matplotlib.pyplot as plt

def sinc(x):
    """Sinc function as used in MATLAB (normalized by pi)."""
    return np.sinc(x / np.pi)

def range_angle(theta):
    """Ensure theta is in the range (-π, π]."""
    while theta > np.pi:
        theta -= 2 * np.pi
    while theta <= -np.pi:
        theta += 2 * np.pi
    return theta

def biarc(x0, y0, th0, x1, y1, th1):
    """
    Compute the biarc passing through (x0, y0) with angle th0
    to (x1, y1) with angle th1.

    Returns:
        l0, theta0, kappa0 - first arc parameters
        l1, theta1, kappa1 - second arc parameters
        xs, ys, thetas - splitting point and intermediate angle
    """
    dx = x1 - x0
    dy = y1 - y0
    d = np.hypot(dx, dy)
    omega = np.arctan2(dy, dx)

    theta0 = omega + range_angle(th0 - omega)
    theta1 = omega + range_angle(th1 - omega)

    dt = (theta1 - theta0) / 2
    t = d * sinc(dt / 2) / sinc(dt)

    thetas = 2 * omega - (theta0 + theta1) / 2
    dt0 = (thetas - theta0) / 2
    dt1 = (thetas - theta1) / 2

    l0 = t / (2 * sinc(dt0))
    l1 = t / (2 * sinc(dt1))
    
    # kappa is the radius of curvature, i.e 1/radius
    kappa0 = 4 * np.sin(dt0) / t
    kappa1 = -4 * np.sin(dt1) / t

    xs = x0 + (t / 2) * np.cos((thetas + theta0) / 2)
    ys = y0 + (t / 2) * np.sin((thetas + theta0) / 2)

    return l0, theta0, kappa0, l1, theta1, kappa1, xs, ys, thetas


if __name__ == "__main__":
    # Example usage:
        
    jCase = 1
    if 1==jCase:    
      x0, y0, th0 = 0, 0, np.pi/4
      outStr='case 1'
    elif 2==jCase:
      x0, y0, th0 = 0, 0, np.pi/8
      outStr ='case 2'
    x1, y1, th1 = 2, 2, np.pi/2
    
    # result = biarc(x0, y0, th0, x1, y1, th1)
    [l0, theta0, kappa0, l1, theta1, kappa1, xs, zs, thetas] =  biarc(x0, y0, th0, x1, y1, th1)
    
    
    print("Biarc parameters:\n")
    print("[l0,        theta0,        kappa0,        l1,        theta1,        kappa1,        xs,        ys,        thetas]")
    outStr=f"{l0:9.3f}, { degrees(theta0):9.3f}, { kappa0:9.3f}, { l1:9.3f}, {degrees(theta1):9.3f}, { kappa1:9.3f}, { xs:9.3f}, { zs:9.3f}, {degrees(thetas):9.3f}"
    print(outStr)
    
    plt.figure()  # Create a new figure
    biarc_plot(x0, y0, l0, theta0, kappa0,
                   x1, y1, l1, theta1, kappa1,
                   fmt1=None, fmt2=None)
    plt.plot(xs, zs, 'go')
    pi2= 0.5*pi
    # pi2 = th0
    # pi2 = 0


    xTxt = xs + cos(thetas + pi2)/kappa0
    zTxt = zs + sin(thetas + pi2)/kappa0
    
    x1, z1 = [xs, xs + cos(thetas + pi2)/kappa0], [zs, zs + sin(thetas + pi2)/kappa0]
    plt.plot(x1, z1)
    x2, z2 = [xs, xs + cos(thetas + pi2)/kappa1], [zs, zs + sin(thetas + pi2)/kappa1]
    plt.plot(x2, z2)    
    x3, z3 = [xs + cos(thetas + pi2)/kappa0, 0], [zs + sin(thetas  + pi2)/kappa0, 0]
    plt.plot(x3, z3)        
    
    plt.text(xTxt - 5, zTxt, outStr)
    plt.show()  
    
    k=0