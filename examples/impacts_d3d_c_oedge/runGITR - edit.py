import os
import subprocess
import shutil
import numpy as np
import libconf
import io
import impacts

import netCDF4 as nc
plasmaProfiles = "sasvw-pa-fav/input/plasmaProfiles.nc"
wallFile = "sasvw-pa-fav/setup/assets/gitr_rz.txt"
ds = nc.Dataset(plasmaProfiles)

#import wall geometry to plot over
with open(wallFile, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])
  
r1 = R[:-1]
r2 = R[1:]
z1 = Z[:-1]
z2 = Z[1:]
slope = np.zeros(len(r1))
Alpha = np.zeros(len(r1))
for i in range(len(r1)):
    if (r2[i]-r1[i])!=0:
        slope[i] = (z2[i]-z1[i])/(r2[i]-r1[i])
        #define angle between material wall and major radius, x
        Alpha[i] = np.abs(np.arctan((z2[i]-z1[i]) / (r2[i]-r1[i])))
    elif (r2[i]-r1[i])==0:
        slope[i] = 100
        Alpha[i] = 89.999*np.pi/180

Alpha = np.abs(Alpha) 

filename = "input/gitrInput.cfg"
angles = ds['Bangle_outer_target'][:]
Bmag = 2.2970
Btheta = 90.5289;
v_parallel = 24561; #24561;
bg = ["C"] #,"He"];
ZZ = [6] #,2];
amu = [12] #,4];
for k in range(0,int(ZZ[0])): #plasma species of bg
    for i in range(0,len(angles)):
    
        with io.open(filename) as f:
            config = libconf.load(f)
            config['backgroundPlasmaProfiles']['Z'] = ZZ[k]
            config['backgroundPlasmaProfiles']['amu'] = amu[k]
            config['backgroundPlasmaProfiles']['Bfield']['r'] = Bmag*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['Bfield']['z'] = Bmag*np.cos(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['Bfield']['y'] = Bmag*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVr'] = -v_parallel*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -v_parallel*np.cos(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVy'] = -v_parallel*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
    
        with io.open(filename, 'w') as f:
            libconf.dump(config, f)
        
        impacts.setup_case(Zi = 6, amu_i = 12, charge_i = 4,
               Zb = ZZ[k], amu_b = amu[k], charge_b = 1,
               Ti = 19.334,
               ni = 4.68e19, Te = 18, ne = 4.68e19, Bmag = 2.297, phi = 90+angles[i],
               nP = 1e3, nT = 1e4, dt = 1.0e-7, height_factor = 4,
               filename = "input/gitrInput.cfg")
        subprocess.run("../build/GITR", shell=True, check=True)
        #shutil.copyfile("output/surface.nc","surface_C"+str(j)+"_loc_"+str(i) +".nc")
        shutil.copyfile("output/positions.nc","positions_"+str(i)+"_" +bg[k] + ".nc")
        os.remove("output/surface.nc")
        os.remove("output/positions.nc")
        os.remove("output/positions.m")
        os.remove("output/particleSource.nc")



















