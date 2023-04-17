# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:58:31 2023

@author: elicr
"""

import numpy as np
import scipy.interpolate as scii
import matplotlib.pyplot as plt
import netCDF4

surface = netCDF4.Dataset("output/surface.nc", "r", format="NETCDF4")
# nP = len(history.dimensions['nP'])
surf = len(surface.dimensions['nSurfaces'])
eng = len(surface.dimensions['nEnergies'])
ang = len(surface.dimensions['nAngles'])
ReflDist = surface.variables['surfReflDist']
SputtDist = surface.variables['surfSputtDist']
EDist = surface.variables['surfEDist']
spylCounts = surface.variables['spylCounts']

x=np.empty([eng,ang])
for i in range(ang):
    x[i] = EDist[0][:][i]


Part = netCDF4.Dataset("output/particleSource.nc", "r", format="NETCDF4")
nP = len(Part.dimensions['nP'])
p = Part.variables['x']
j = Part.variables['y']
k = Part.variables['z']

x = np.empty(nP)
y = np.empty(nP)
z = np.empty(nP)

for i in range(nP):
    x[i] = p[i]
    y[i] = j[i]
    z[i] = k[i]

# Writing geometry to txt
n = int(len(x))
coord = np.empty([n,2])
for i in range(n):
    coord[i,:] = [round(x[i],6),round(z[i],6)]    
myFile = open('exp_rz.txt', 'r+')
np.savetxt(myFile, coord)
myFile.close()