#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import netCDF4
import numpy as np
import scipy.interpolate as scii
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

case = 1
device = 'local'

plt.rcParams.update({'font.size':12})
plt.rcParams.update({'lines.linewidth':2})

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
    W_indices = np.arange(11,22)
elif case == 2: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-pa-unfav/setup'
    W_indices = np.arange(11,22)
elif case == 3: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-vertex-fav/setup'
    W_indices = np.arange(16,25)
elif case == 4: 
    setup_directory = host_dir+'/GITR_processing/examples/sasvw-vertex-unfav/setup'
    W_indices = np.arange(16,25)
else:
    print('Error: case number must be 1, 2, 3, or 4')

#import refined rmrs at the W surface
rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt'

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

species_list = np.array(['D1','C1','C2','C3','C4','C5','C6'])
angles = np.arange(0,90,1)

##############################################################################
# animate flux as a function of energy across the W tile and
# plot average energy across the W tile
##############################################################################

def calculate_flux_over_energy_AND_avg_energy(verbose=1):
    energies = np.zeros((len(rmrsFine), len(species_list), 240))
    
    total_relative_flux = np.zeros((len(rmrsFine), len(species_list)))
    relative_flux_over_energy = np.zeros((len(rmrsFine), len(species_list), 240))
    normalized_flux_over_energy = np.zeros((len(rmrsFine), len(species_list), 240))
    
    average_energy = np.zeros((len(species_list), len(rmrsFine)))
    
    for i in range(len(rmrsFine)):
        if verbose: print('\nLine Segment:', i, 'of', len(rmrsFine)-1)
        pathname = host_dir+'/hpic2-data/case'+str(case)+'/IEAD_data/x'+str(i+1)
        
        for s in range(len(species_list)):
            species = species_list[s]
            if verbose: print('Species:', species)
            
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
            
            energies[i,s] = np.linspace(dE, Emax, 240)
            
            ##########################################################################
            # Calculate total relative flux for normalization
            ##########################################################################
            
            # all relative flux data is in units of # of particles over energy and angle grids
            dist = np.loadtxt(filename, dtype=float)
            
            if np.sum(dist) != 0:
                
                # calculate relative flux for each energy bin, normalized by the total relative flux
                relative_flux_over_energy[i,s] = np.sum(dist, axis=1)
                total_relative_flux[i,s] = np.sum(dist)
                normalized_flux_over_energy[i,s] = relative_flux_over_energy[i,s]/total_relative_flux[i,s]
                
                # use normalized fluxes to find the weighted average energy
                average_energy[s,i] = np.sum(normalized_flux_over_energy[i,s] * energies[i,s])
                
    return energies, normalized_flux_over_energy, average_energy

def plot_avg_energy():
    # get estimated energy using Stangeby method to plot against
    profilesFile = setup_directory+'/../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices]
    teCoarse = profiles.variables['te_inner_target'][W_indices]
    tiCoarse = profiles.variables['te_inner_target'][W_indices]
    fte = scii.interp1d(rmrsCoarse, teCoarse, fill_value='extrapolate')
    fti = scii.interp1d(rmrsCoarse, tiCoarse, fill_value='extrapolate')
    teFine = fte(rmrsFine)
    tiFine = fti(rmrsFine)
    estimated_energy = np.zeros((len(species_list), len(rmrsFine)))
    estimated_energy[0] = 3*(1)*teFine + 2*tiFine
    for s in range(1,len(species_list)):
        estimated_energy[s] = 3*(s+1)*teFine + 2*tiFine
    
    # get hPIC2-calculated weighted average incident energies
    energies, normalized_flux_over_energy, average_energy = calculate_flux_over_energy_AND_avg_energy(verbose=0)
    
    plt.close()
    plt.hlines(45, np.min(rmrsFine), np.max(rmrsFine), color='gray')
    colors = ['black', 'firebrick', 'darkorange', 'gold', 'limegreen', 'dodgerblue', 'mediumpurple']
    for s in range(len(species_list)):
        plt.plot(rmrsFine, average_energy[s], color=colors[s], label=species_list[s])
        plt.plot(rmrsFine, estimated_energy[s], '--', color=colors[s])
    
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Energy [eV]')
    plt.title('Average Incident Energy\n Calculated by hPIC2')
    if case==4: 
        plt.legend(loc=4)
    else: 
        plt.legend()
    #plt.yscale('log')
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(1,3))
    plt.savefig(setup_directory+'/plots/surface-profiles/hPIC2_avg_incident_energy.png')
    
if __name__ == "__main__":
    #plot_avg_energy()
    
    energies, normalized_flux_over_energy, average_energy = calculate_flux_over_energy_AND_avg_energy(verbose=0)
    ymin = np.min(normalized_flux_over_energy)
    ymax = np.max(normalized_flux_over_energy)
    plt.close()
    fig, ax = plt.subplots()
    colors = ['black', 'firebrick', 'darkorange', 'gold', 'limegreen', 'dodgerblue', 'mediumpurple']
    
    def update(i):
        ax.cla()
        ax.set_xlim([np.min(energies), np.max(energies)])
        ax.set_ylim([ymin, ymax])
        plt.xlabel('Energy [eV]')
        plt.ylabel('Normalized relative flux [arb]')
        plt.title('D-Dsep = '+str(rmrsFine[i]))
        for s in range(len(species_list)):
            ax.plot(energies[i,s], normalized_flux_over_energy[i,s], color=colors[s])
            plt.vlines(average_energy[s,i], ymin, ymax, color=colors[s])
        return
    
    ani = FuncAnimation(fig, update, frames=len(rmrsFine), blit=False)
    ani.save(filename=setup_directory+'/plots/particle-source/hPIC2_energy_distribution_animated.gif',\
             writer="pillow")
    