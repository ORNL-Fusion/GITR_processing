#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import netCDF4
import numpy as np
import scipy.interpolate as scii
import pandas as pd

case = 1
species = 'C1'
device = 'perlmutter'

##############################################################################
# path and file definitions
##############################################################################

if device == 'local':
    host_dir = '/Users/Alyssa/Dev'
elif device == 'perlmutter':
    host_dir = '/pscratch/sd/h/hayes'
else: 
    print('Error: setup_directories must be defined for this device')

if case == 1: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-pa-fav/setup'
elif case == 2: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-pa-unfav/setup'
elif case == 3: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-vertex-fav/setup'
elif case == 4: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-vertex-unfav/setup'
else:
    print('Error: case number must be 1, 2, 3, or 4')

#file definitions
rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt'
if 'C' in species:
    ftBfile = setup_directory+'/assets/ftridynBackgroundC.nc'
elif 'D' in species:
    ftBfile = setup_directory+'/assets/ftridynBackgroundD.nc'
else:
    print('Error: species must contain D or C')

#import refined rmrs at the W surface
with open(rmrs_fine_file, 'r') as file:
    rmrs_fine = file.readlines()   
rmrsFine = np.array(rmrs_fine,dtype='float')

def complete_filename(species, path):
    result = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if species in f:
                result=f
    return path+'/'+result

def get_ft_spyld(S, surfE, surfA, ftBfile):
    #import sputtering yield tables for incident ions on W
    ftB = netCDF4.Dataset(ftBfile, "r", format="NETCDF4")
    
    try:
        spyld = ftB.variables['spyld'][S][:]
        ftE = ftB.variables['E'][:]
        ftA = ftB.variables['A'][:]
        
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA), bounds_error=False, fill_value=None)
        
    except:
        spyld = ftB.variables['spyld'][S][:]
        ftE = ftB.variables['E'][:]
        ftA = ftB.variables['A'][:]
        
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA), bounds_error=False, fill_value=None)
    
    return surfY

##############################################################################
# MAIN EXECUTED CODE
##############################################################################

angles = np.arange(0,90,1)

total_relative_flux = np.zeros(len(rmrsFine))
relative_spyld = np.zeros(len(rmrsFine))
effective_spyld = np.zeros(len(rmrsFine))

for i in range(len(rmrsFine)):
    print('\nLine Segment:', i, 'of', len(rmrsFine)-1)
    pathname = host_dir+'/hpic2-data/case1/IEAD_data/x'+str(i+1)
    
    # autofill filenames for each species
    filename = complete_filename('_'+species+'_', pathname)
    
    ##########################################################################
    # Calculate energy grid for each species because Emax is species-dependent
    ##########################################################################
    
    # pull dE from the filename
    Emarker = filename.index('_dEeV_')
    
    dE = float(filename[Emarker+6:Emarker+13])
    
    # calculate Emax and the energy spectrum from the filename because there are always 240 bins
    Emax = 240 * dE
    
    energies = np.arange(dE, Emax, dE)
    
    ##########################################################################
    # Interpolate the relative sputtering yield from the relative fluxes for
    # all species at each line segment by loop through all (E,A) pairs
    ##########################################################################
    
    # all relative flux data is in units of # of particles over energy and angle grids
    dist = np.loadtxt(filename, dtype=float)
    
    # calculate total relative flux for later normalization
    if np.sum(dist) != 0:
        total_relative_flux[i] = np.sum(dist)
        
        # calculate relative sputtering yields
        for A in range(len(angles)): 
            print('Angle:', angles[A], 'of', angles[-1])
            
            for E in range(len(energies)): 
                relative_flux = dist[E, A]
                
                spyld = get_ft_spyld(0, energies[E], angles[A], ftBfile)
                
                relative_spyld[i] +=  relative_flux * spyld


    # normalize relative sputtering yields by dividing by the total relative flux
    effective_spyld[i] = relative_spyld[i] / total_relative_flux[i]
    print('Effective spyld:', effective_spyld[i])

##########################################################################
# save effective sputtering yields to a file
##########################################################################

eff_spyld_df = pd.DataFrame(np.transpose(effective_spyld), \
                            columns = [species])

eff_spyld_df.to_csv(setup_directory+'/assets/eff_spyld_hPIC_case'+str(case)+'_'+species+'.csv', index=False)

