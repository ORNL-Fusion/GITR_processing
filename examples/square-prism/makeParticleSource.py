#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import shutil
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

def line_source(nP = int(2e2), plotting=1):
    x1 = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    z1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0])
    x2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 0.0])
    z2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.0])
    slope = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e12])
    inDir = np.array([-1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1])
    plt.close()
    
    for i in range(0,len(x1)):
        if slope[i]==0: 
            perpSlope = 1.0e12
        else:
            perpSlope = -np.sign(slope[i])/np.abs(slope[i])

        rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1)
        zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp)
        if plotting: plt.quiver([x1[i] + (x2[i]-x1[i])/2], [z1[i] + (z2[i]-z1[i])/2], [rPerp/10], [zPerp/10], width=0.0015, scale=5, headwidth=4)
    
    x_23 = 2.5*np.ones(int(nP/4))
    x_34 = 3.5*np.ones(int(nP/4))
    x_45 = 4.5*np.ones(int(nP/4))
    x_56 = 5.5*np.ones(int(nP/4))
    
    x = np.hstack((x_23, x_34, x_45, x_56))
    y = np.zeros(nP)
    z = (1e-7)*np.ones(nP)
    
    conversion = np.sqrt(2 * 1.6021773e-19 / 1.6605E-27 / 183) #eV to m/s
    vx = np.zeros(nP)*conversion
    vy = np.zeros(nP)*conversion
    vz = 200*np.ones(nP)*conversion

    #########################################
    #make NetCDF Particle Source file
    #########################################

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nP)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = x
    yyy[:] = y
    zzz[:] = z
    vxx[:] = vx
    vyy[:] = vy
    vzz[:] = vz
    rootgrp.close()
    shutil.move('particleSource.nc', 'input/particleSource.nc')



if __name__ == "__main__":
    line_source()