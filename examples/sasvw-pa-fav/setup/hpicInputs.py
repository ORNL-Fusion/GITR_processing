###############################################################################
# READ ME
# 
# It is best to use an .eq or an .equ file for the Bfield whenever possible
# For this particular case, the Bfield components are relatively constant
###############################################################################

import numpy as np
import scipy.interpolate as scii
import netCDF4

###############################################################################
# Import coarse variable dataset from processed SOLPS-ITER output
###############################################################################

profiles_file = '../input/plasmaProfiles.nc'
profiles = netCDF4.Dataset(profiles_file)
W_indices_coarse = np.arange(11,22)

###############################################################################
# Import refined points of interest
###############################################################################

rmrs_fine_file = '../setup/assets/rmrs_fine.txt'

with open(rmrs_fine_file, 'r') as file:
    rmrs_fine = file.readlines()   
rmrsFine = np.array(rmrs_fine,dtype='float')

###############################################################################
# Parse coarse dataset for desired variables: Bmag, Bangle, te, ti, ne, ni(s)
###############################################################################

rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices_coarse]
BmagCoarse = profiles.variables['Bmag_inner_target'][W_indices_coarse]
BangleCoarse = 180-profiles.variables['Bangle_inner_target'][W_indices_coarse] # w.r.t. normal incidence
teCoarse = profiles.variables['te_inner_target'][W_indices_coarse]
tiCoarse = profiles.variables['ti_inner_target'][W_indices_coarse]
neCoarse = profiles.variables['ne_inner_target'][W_indices_coarse]

niCoarseD0 = np.abs(profiles.variables['flux_inner_target'][0][W_indices_coarse])
niCoarseD = np.abs(profiles.variables['flux_inner_target'][1][W_indices_coarse])
niCoarseC0 = np.abs(profiles.variables['flux_inner_target'][2][W_indices_coarse])
niCoarseC1 = np.abs(profiles.variables['flux_inner_target'][3][W_indices_coarse])
niCoarseC2 = np.abs(profiles.variables['flux_inner_target'][4][W_indices_coarse])
niCoarseC3 = np.abs(profiles.variables['flux_inner_target'][5][W_indices_coarse])
niCoarseC4 = np.abs(profiles.variables['flux_inner_target'][6][W_indices_coarse])
niCoarseC5 = np.abs(profiles.variables['flux_inner_target'][7][W_indices_coarse])
niCoarseC6 = np.abs(profiles.variables['flux_inner_target'][8][W_indices_coarse])

###############################################################################
# Define linear interpolation functions
###############################################################################

BmagFunc = scii.interp1d(rmrsCoarse, BmagCoarse, fill_value='extrapolate')
BangleFunc = scii.interp1d(rmrsCoarse, BangleCoarse, fill_value='extrapolate')
teFunc = scii.interp1d(rmrsCoarse, teCoarse, fill_value='extrapolate')
tiFunc = scii.interp1d(rmrsCoarse, tiCoarse, fill_value='extrapolate')
neFunc = scii.interp1d(rmrsCoarse, neCoarse, fill_value='extrapolate')

niFuncD0 = scii.interp1d(rmrsCoarse,niCoarseD0,fill_value='extrapolate')
niFuncD = scii.interp1d(rmrsCoarse,niCoarseD,fill_value='extrapolate')
niFuncC0 = scii.interp1d(rmrsCoarse,niCoarseC0,fill_value='extrapolate')
niFuncC1 = scii.interp1d(rmrsCoarse,niCoarseC1,fill_value='extrapolate')
niFuncC2 = scii.interp1d(rmrsCoarse,niCoarseC2,fill_value='extrapolate')
niFuncC3 = scii.interp1d(rmrsCoarse,niCoarseC3,fill_value='extrapolate')
niFuncC4 = scii.interp1d(rmrsCoarse,niCoarseC4,fill_value='extrapolate')
niFuncC5 = scii.interp1d(rmrsCoarse,niCoarseC5,fill_value='extrapolate')
niFuncC6 = scii.interp1d(rmrsCoarse,niCoarseC6,fill_value='extrapolate')

###############################################################################
# Plug refined rmrs points into interpolation functions
###############################################################################

BmagFine = BmagFunc(rmrsFine)
BangleFine = BangleFunc(rmrsFine)
teFine = teFunc(rmrsFine)
tiFine = tiFunc(rmrsFine)
neFine = neFunc(rmrsFine)

niFineD0 = niFuncD0(rmrsFine)
niFineD = niFuncD(rmrsFine)
niFineC0 = niFuncC0(rmrsFine)
niFineC1 = niFuncC1(rmrsFine)
niFineC2 = niFuncC2(rmrsFine)
niFineC3 = niFuncC3(rmrsFine)
niFineC4 = niFuncC4(rmrsFine)
niFineC5 = niFuncC5(rmrsFine)
niFineC6 = niFuncC6(rmrsFine)