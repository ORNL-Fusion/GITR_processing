#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import netCDF4
import numpy as np
import scipy.interpolate as scii
import pandas as pd

case = 1
if case == 0: case = int(sys.argv[1])

device = 'local'

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
    print('Error: Case number must be 1, 2, 3, or 4')

#file definitions
rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt'
ftDfile = setup_directory+'/assets/ftridynBackgroundD.nc'
ftCfile = setup_directory+'/assets/ftridynBackgroundC.nc'

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
        
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA), bounds_error=False)
        
    except:
        spyld = ftB.variables['spyld'][S][:]
        ftE = ftB.variables['E'][:]
        ftA = ftB.variables['A'][:]
        
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA), bounds_error=False)
    
    return surfY

##############################################################################
# MAIN EXECUTED CODE
##############################################################################

angles = np.arange(0,90,1)

total_relative_fluxD = np.zeros(len(rmrsFine))
total_relative_fluxC1 = np.zeros(len(rmrsFine))
total_relative_fluxC2 = np.zeros(len(rmrsFine))
total_relative_fluxC3 = np.zeros(len(rmrsFine))
total_relative_fluxC4 = np.zeros(len(rmrsFine))
total_relative_fluxC5 = np.zeros(len(rmrsFine))
total_relative_fluxC6 = np.zeros(len(rmrsFine))

relative_spyldD = np.zeros(len(rmrsFine))
relative_spyldC1 = np.zeros(len(rmrsFine))
relative_spyldC2 = np.zeros(len(rmrsFine))
relative_spyldC3 = np.zeros(len(rmrsFine))
relative_spyldC4 = np.zeros(len(rmrsFine))
relative_spyldC5 = np.zeros(len(rmrsFine))
relative_spyldC6 = np.zeros(len(rmrsFine))

for i in range(3,len(rmrsFine)):
    print('\nLine Segment:', i, 'of', len(rmrsFine))
    pathname = host_dir+'/hpic2-data/case1/IEAD_data/x'+str(i+1)
    
    # autofill filenames for each species
    filenameD = complete_filename('_D1_', pathname)
    filenameC1 = complete_filename('_C1_', pathname)
    filenameC2 = complete_filename('_C2_', pathname)
    filenameC3 = complete_filename('_C3_', pathname)
    filenameC4 = complete_filename('_C4_', pathname)
    filenameC5 = complete_filename('_C5_', pathname)
    filenameC6 = complete_filename('_C6_', pathname)
    
    ##########################################################################
    # Calculate energy grid for each species because Emax is species-dependent
    ##########################################################################
    
    # pull dE from the filename
    Emarker_D = filenameD.index('_dEeV_')
    Emarker_C1 = filenameC1.index('_dEeV_')
    Emarker_C2 = filenameC2.index('_dEeV_')
    Emarker_C3 = filenameC3.index('_dEeV_')
    Emarker_C4 = filenameC4.index('_dEeV_')
    Emarker_C5 = filenameC5.index('_dEeV_')
    Emarker_C6 = filenameC6.index('_dEeV_')
    
    dE_D = float(filenameD[Emarker_D+6:Emarker_D+13])
    dE_C1 = float(filenameC1[Emarker_C1+6:Emarker_C1+13])
    dE_C2 = float(filenameC2[Emarker_C2+6:Emarker_C2+13])
    dE_C3 = float(filenameC3[Emarker_C3+6:Emarker_C3+13])
    dE_C4 = float(filenameC4[Emarker_C4+6:Emarker_C4+13])
    dE_C5 = float(filenameC5[Emarker_C5+6:Emarker_C5+13])
    dE_C6 = float(filenameC6[Emarker_C6+6:Emarker_C6+13])
    
    # calculate Emax and the energy spectrum from the filename because there are always 240 bins
    Emax_D = 240 * dE_D
    Emax_C1 = 240 * dE_C1
    Emax_C2 = 240 * dE_C2
    Emax_C3 = 240 * dE_C3
    Emax_C4 = 240 * dE_C4
    Emax_C5 = 240 * dE_C5
    Emax_C6 = 240 * dE_C6
    
    energiesD = np.arange(dE_D, Emax_D, dE_D)
    energiesC1 = np.arange(dE_C1, Emax_C1, dE_C1)
    energiesC2 = np.arange(dE_C2, Emax_C2, dE_C2)
    energiesC3 = np.arange(dE_C3, Emax_C3, dE_C3)
    energiesC4 = np.arange(dE_C4, Emax_C4, dE_C4)
    energiesC5 = np.arange(dE_C5, Emax_C5, dE_C5)
    energiesC6 = np.arange(dE_C6, Emax_C6, dE_C6)
    
    '''
    # set E=1e3 eV for all E>1e3 eV
    # we anticipate a factor of 2 at max between the peak in the spyld curve and E=1e3 eV
    # citation: Eckstein book, "Sputtering by Particle Bombardment"
    if np.max(energiesD)>1000:
        energiesD[np.where(energiesD>1000)[0]] = 1000
    if np.max(energiesC1)>1000:
        energiesC1[np.where(energiesC1>1000)[0]] = 1000
    if np.max(energiesC2)>1000:
        energiesC2[np.where(energiesC2>1000)[0]] = 1000
    if np.max(energiesC3)>1000:
        energiesC3[np.where(energiesC3>1000)[0]] = 1000
    if np.max(energiesC4)>1000:
        energiesC4[np.where(energiesC4>1000)[0]] = 1000
    if np.max(energiesC5)>1000:
        energiesC5[np.where(energiesC5>1000)[0]] = 1000
    if np.max(energiesC6)>1000:
        energiesC6[np.where(energiesC6>1000)[0]] = 1000
    '''
    
    ##########################################################################
    # Interpolate the relative sputtering yield from the relative fluxes for
    # all species at each line segment by loop through all (E,A) pairs
    ##########################################################################
    
    # all relative flux data is in units of # of particles over energy and angle grids
    distD = np.loadtxt(filenameD, dtype=float)
    distC1 = np.loadtxt(filenameC1, dtype=float)
    distC2 = np.loadtxt(filenameC2, dtype=float)
    distC3 = np.loadtxt(filenameC3, dtype=float)
    distC4 = np.loadtxt(filenameC4, dtype=float)
    distC5 = np.loadtxt(filenameC5, dtype=float)
    distC6 = np.loadtxt(filenameC6, dtype=float)
    
    # calculate total relative flux for later normalization
    total_relative_fluxD[i] = np.sum(distD)
    total_relative_fluxC1[i] = np.sum(distC1)
    total_relative_fluxC2[i] = np.sum(distC2)
    total_relative_fluxC3[i] = np.sum(distC3)
    total_relative_fluxC4[i] = np.sum(distC4)
    total_relative_fluxC5[i] = np.sum(distC5)
    total_relative_fluxC6[i] = np.sum(distC6)
    
    # calculate relative sputtering yields
    for A in range(len(angles)): 
        print('Angle:', angles[A], 'of', angles[-1])
        
        for E in range(len(energiesD)): 
            relative_fluxD = distD[E, A]
            relative_fluxC1 = distC1[E, A]
            relative_fluxC2 = distC2[E, A]
            relative_fluxC3 = distC3[E, A]
            relative_fluxC4 = distC4[E, A]
            relative_fluxC5 = distC5[E, A]
            relative_fluxC6 = distC6[E, A]
            
            spyldD = get_ft_spyld(1, energiesD[E], angles[A], ftDfile)
            spyldC1 = get_ft_spyld(0, energiesC1[E], angles[A], ftCfile)
            spyldC2 = get_ft_spyld(0, energiesC2[E], angles[A], ftCfile)
            spyldC3 = get_ft_spyld(0, energiesC3[E], angles[A], ftCfile)
            spyldC4 = get_ft_spyld(0, energiesC4[E], angles[A], ftCfile)
            spyldC5 = get_ft_spyld(0, energiesC5[E], angles[A], ftCfile)
            spyldC6 = get_ft_spyld(0, energiesC6[E], angles[A], ftCfile)
        
            relative_spyldD[i] +=  relative_fluxD * spyldD
            relative_spyldC1[i] +=  relative_fluxC1 * spyldC1
            relative_spyldC2[i] +=  relative_fluxC2 * spyldC2
            relative_spyldC3[i] +=  relative_fluxC3 * spyldC3
            relative_spyldC4[i] +=  relative_fluxC4 * spyldC4
            relative_spyldC5[i] +=  relative_fluxC5 * spyldC5
            relative_spyldC6[i] +=  relative_fluxC6 * spyldC6


# normalize relative sputtering yields by dividing by the total relative flux

effective_spyldD = relative_spyldD / total_relative_fluxD
effective_spyldC1 = relative_spyldC1 / total_relative_fluxC1
effective_spyldC2 = relative_spyldC2 / total_relative_fluxC2
effective_spyldC3 = relative_spyldC3 / total_relative_fluxC3
effective_spyldC4 = relative_spyldC4 / total_relative_fluxC4
effective_spyldC5 = relative_spyldC5 / total_relative_fluxC5
effective_spyldC6 = relative_spyldC6 / total_relative_fluxC6

##########################################################################
# save effective sputtering yields to a file
##########################################################################

eff_spyld_df = pd.DataFrame(np.transpose(np.vstack((effective_spyldD, effective_spyldC1, effective_spyldC2, \
                                  effective_spyldC3, effective_spyldC4, effective_spyldC5, effective_spyldC6))), \
                            columns = ['D','C1','C2','C3','C4','C5','C6'])

eff_spyld_df.to_csv(setup_directory+'/assets/eff_spyld_hPIC_case'+case+'.csv', index=False)

