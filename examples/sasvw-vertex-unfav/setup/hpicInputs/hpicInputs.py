###############################################################################
# READ ME
# Author: Alyssa Hayes
# 
# This script was originally intended for use by Andy Liu and Josh Hoffman. 
# This script returns variables from GITR pre-processing for use as hPIC2 inputs. 
# All plasma parameters Bmag, te, ti, ne, and ni are sheath entrance parameters. 
# Each 'ith' variable set corresponds to one refined line segment on a W surface. 
# 
# Variables returned from 2D bfield.nc: Bmag and Bangle
# Variables returned from 2D plasmaProfiles.nc: ne, te, ti
# Variables returned from 1D plasmaProfiles.nc: ni(ns)
###############################################################################

import numpy as np
import scipy.interpolate as scii
import netCDF4
import matplotlib.pyplot as plt

case_number = 4
W_indices_coarse = np.arange(16,25)

rmrs_fine_file = 'rmrs_fine.txt'
gitr_rz = 'gitr_rz.txt'
W_fine_file = 'W_fine.txt'

profiles_file = 'plasmaProfiles.nc'
profiles = netCDF4.Dataset(profiles_file)
bfield_file = 'bField.nc'
bfield = netCDF4.Dataset(bfield_file)

# toggle this to see test comparisons of input vars
plot_blocker = True

###############################################################################
# Import refined rmrs points of interest for interpolation off inner_target
###############################################################################

with open(rmrs_fine_file, 'r') as file:
    rmrs_fine = file.readlines()   

rmrsFine = np.array(rmrs_fine,dtype='float')
print('\nThese should all match')
print('Length of W surface by rmrs:', len(rmrsFine))

###############################################################################
# Import refined (r,z) points of interest for interpolation off full 2D grid
###############################################################################

# Import full wall (R,Z)
with open(gitr_rz, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])
    
# Slice out only the W indices to get only the wall (R,Z) for W surfaces
with open(W_fine_file, 'r') as file:
    W_fine = file.readlines()
W_fine = np.array(W_fine,dtype='int')

R = R[W_fine]
Z = Z[W_fine]

# Calculate (R,Z) midpoints and check lengths
R1 = R[:-1]
R2 = R[1:]
Z1 = Z[:-1]
Z2 = Z[1:]
dR = R2-R1
dZ = Z2-Z1
R_midpoints_fine = (R1+R2)/2
Z_midpoints_fine = (Z1+Z2)/2
print('Length of W surface by (R,Z):', len(R_midpoints_fine), len(Z_midpoints_fine),'\n')    

###############################################################################
# Get (Br, Bt, Bz) from bField.nc 2D grid interpolation
###############################################################################

gridz = bfield.variables['z'][:]
gridr = bfield.variables['r'][:]

br_2D = bfield.variables['br'][:]
bt_2D = bfield.variables['bt'][:]
bz_2D = bfield.variables['bz'][:]

br = scii.interpn((gridz,gridr),br_2D,(Z_midpoints_fine,R_midpoints_fine))
bt = scii.interpn((gridz,gridr),bt_2D,(Z_midpoints_fine,R_midpoints_fine))
bz = scii.interpn((gridz,gridr),bz_2D,(Z_midpoints_fine,R_midpoints_fine))

###############################################################################
# Get Bmag from (br, bt, bz)
#
# Get Bangle using the law that cos(θ) = n•B/(|n||B|)
# Where n is the vector representing the normal direction from the surface
# such that n = -dZ i + 0 j + dR k
###############################################################################

top = -dZ*br + dR*bz
Bmag = np.sqrt(br**2 + bt**2 + bz**2)
nmag = np.sqrt(dZ**2 + dR**2)
bottom = Bmag * nmag
Bangle = 180-np.rad2deg(np.arccos(top/bottom)) # w.r.t. normal incidence

# Test Bangle comparison against the Bangle from the 1D plasmaProfiles.nc var
rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices_coarse]
Bangle_coarse = 180-profiles.variables['Bangle_inner_target'][W_indices_coarse] # w.r.t. normal incidence
Bangle_func = scii.interp1d(rmrsCoarse, Bangle_coarse, fill_value='extrapolate')
Bangle_test = Bangle_func(rmrsFine)

plt.close()
plt.plot(rmrsFine, Bmag, color='black')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Bfield Strength [T]')
plt.title('Case '+str(case_number)+' Bmag')
plt.show(block=plot_blocker)

plt.close()
plt.plot(rmrsFine, Bangle_test, color='pink', label='from "inner_target"')
plt.plot(rmrsFine, Bangle, color='black', label='interpolated from 2D bfield.nc')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Degrees from normal incidence')
plt.title('Case '+str(case_number)+' Bangle')
plt.legend()
plt.show(block=plot_blocker)

###############################################################################
# Get ne, te, ti from profiles.nc auto-interpolated vars_inner_target
###############################################################################

rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices_coarse]
ne_coarse = profiles.variables['ne_inner_target'][W_indices_coarse]
te_coarse = profiles.variables['te_inner_target'][W_indices_coarse]
ti_coarse = profiles.variables['ti_inner_target'][W_indices_coarse]

ne_func = scii.interp1d(rmrsCoarse, ne_coarse, fill_value='extrapolate')
te_func = scii.interp1d(rmrsCoarse, te_coarse, fill_value=te_coarse[0], bounds_error=False)
ti_func = scii.interp1d(rmrsCoarse, ti_coarse, fill_value='extrapolate')

ne = ne_func(rmrsFine)
te = te_func(rmrsFine)
ti = ti_func(rmrsFine)

plt.close()
plt.plot(rmrsFine, ne, color='pink')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Density [m$^{-3}$]')
plt.title('Case '+str(case_number)+' ne')
plt.legend()
plt.show(block=plot_blocker)

plt.close()
plt.plot(rmrsFine, te, color='pink')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Temperature [eV]')
plt.title('Case '+str(case_number)+' te')
plt.show(block=plot_blocker)

plt.close()
plt.plot(rmrsFine, ti, color='pink')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Temperature [eV]')
plt.title('Case '+str(case_number)+' ti')
plt.show(block=plot_blocker)

###############################################################################
# Get ni(ns) from GITR inner target profiles variable
###############################################################################

rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices_coarse]
def extract_ni_by_species(ns):
    ni_coarse = profiles.variables['ni_inner_target'][ns][W_indices_coarse]
    ni_func = scii.interp1d(rmrsCoarse, ni_coarse, fill_value='extrapolate')
    ni = ni_func(rmrsFine)
    return ni

niD0 = extract_ni_by_species(0)
niD1 = extract_ni_by_species(1)
niC0 = extract_ni_by_species(2)
niC1 = extract_ni_by_species(3)
niC2 = extract_ni_by_species(4)
niC3 = extract_ni_by_species(5)
niC4 = extract_ni_by_species(6)
niC5 = extract_ni_by_species(7)
niC6 = extract_ni_by_species(8)

plt.close()
plt.plot(rmrsFine, niD0, color='lightgray', label='D0', linewidth=5) # all 1e10
plt.plot(rmrsFine, niD1, color='rosybrown', label='D1')
plt.plot(rmrsFine, niC0, color='black', label='C0') # all 1e10
plt.plot(rmrsFine, niC1, color='firebrick', label='C1')
plt.plot(rmrsFine, niC2, color='darkorange', label='C2')
plt.plot(rmrsFine, niC3, color='gold', label='C3')
plt.plot(rmrsFine, niC4, color='limegreen', label='C4')
plt.plot(rmrsFine, niC5, color='dodgerblue', label='C5')
plt.plot(rmrsFine, niC6, color='mediumpurple', label='C6')
plt.yscale('log')
plt.xlabel('D-Dsep [m]')
plt.ylabel('Density [m$^{-3}$]')
plt.title('Case '+str(case_number)+' ni')
plt.legend()
plt.show(block=plot_blocker)

###############################################################################
# Returned variable array names as a function of location are below.
# Each index of the array corresponds to a different spot on the surface
# such that each 'ith' index corresponds to the 'ith' hPIC2 simulation. 
#
# For Bmag, te, ti, ne, and ni: assume that these are the values at the
# SHEATH ENTRANCE, not at the surface. SOLPS-ITER uses the Bohm Criterion
# to solve for plasma parameters at the entrance, and does not resolve
# sheath physics into the wall. 
# 
# Use these variable names as written exactly:
# 
# Bmag, Bangle, ne, te, ti, 
# niD0, niD1, niC0, niC1, niC2, niC3, niC4, niC5, niC6
###############################################################################

###############################################################################
# Space for UIUC to write I/O from Python arrays to TOML
###############################################################################










