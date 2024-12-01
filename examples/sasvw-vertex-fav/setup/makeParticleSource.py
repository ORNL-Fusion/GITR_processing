import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import numpy as np
import scipy.interpolate as scii
import matplotlib.pyplot as plt
import netCDF4

import gitr
import solps
import Particles

def init():
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':1.2})
    plt.rcParams.update({'lines.markersize':1})

def point_source(nP = int(2e2)):
    conversion = np.sqrt(2 * 1.6021773e-19 / 1.6605E-27 / 183) #eV to m/s
    
    x = 1.498*np.ones(nP)
    y = np.zeros(nP)
    z = 1.1977*np.ones(nP)
    vx = -0.447*np.ones(nP)*conversion
    vy = 3.5072*np.ones(nP)*conversion
    vz = 3.5355*np.ones(nP)*conversion
    #5 eV, 1 degree off the surface, 83.7363 degree between surface and lab frame

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


def distributed_source(nP, surfW, tile_shift_indices=[], Bangle_shift_indices=[], \
            geom = '../input/gitrGeometry.cfg', \
            profiles_file = '../input/plasmaProfiles.nc', \
            gitr_rz = '../setup/assets/gitr_rz.txt', \
            rmrs_fine_file = '../setup/assets/rmrs_fine.txt', \
            W_fine_file = '../setup/assets/W_fine.txt', \
            ftDFile = '../setup/assets/ftridynBackgroundD.nc', \
            ftCFile = '../setup/assets/ftridynBackgroundC.nc', \
            ftWFile = '../input/ftridynSelf.nc', \
            configuration = 'random', \
            use_fractal_tridyn_outgoing_IEADS = 0, \
            plot_variables = 0):
    
    #import wall geometry to plot over
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()

    #import coarse rmrs at W surface
    profiles = netCDF4.Dataset(profiles_file)
    r_right_target = profiles.variables['r_inner_target']
    z_right_target = profiles.variables['z_inner_target']
    rmrsCoords = profiles.variables['rmrs_inner_target'][surfW]
    rmrsCoarse = profiles.variables['rmrs_inner_target_midpoints'][surfW]

    #import refined rmrs at the W surface
    with open(rmrs_fine_file, 'r') as file:
        rmrs_fine = file.readlines()   
    rmrsFine = np.array(rmrs_fine,dtype='float')
        
    #import W surface indices
    with open(W_fine_file, 'r') as file:
        W_fine = file.readlines()
    W_fine = np.array(W_fine,dtype='int')

    R = np.zeros(len(wall))
    Z = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        R[i] = float(point[0])
        Z[i] = float(point[1])
    
    R = R[W_fine]
    Z = Z[W_fine]
    
    r1 = R[:-1]
    r2 = R[1:]
    z1 = Z[:-1]
    z2 = Z[1:]
    
    if plot_variables==1:
        plt.close()
        plt.plot(r_right_target, z_right_target, '-k', label='Carbon', linewidth=0.5, zorder=0)
       # plt.plot(r_right_target[surfW], z_right_target[surfW], 'purple', label='Profiles Tungsten', linewidth=3, zorder=1)
        plt.plot(R, Z, 'violet', label='Tungsten', linewidth=0.6, zorder=2)
        plt.scatter(R, Z, marker='_', color='violet', s=8, zorder=3)
        plt.axis('scaled')
        plt.legend()
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('Upper Outer SAS-VW Divertor in DIII-D \n makeParticleSource')
        plt.savefig('plots/geom/makePSGeom.png')

    slope = np.zeros(len(r1))
    Alpha = np.zeros(len(r1))
    for i in range(len(r1)):
        if (r2[i]-r1[i])!=0:
            slope[i] = (z2[i]-z1[i])/(r2[i]-r1[i])
            #define angle between material wall and major radius, x
            Alpha[i] = np.abs(np.arctan((z2[i]-z1[i]) / (r2[i]-r1[i])))
        elif (r2[i]-r1[i])==0:
            slope[i] = 100
            Alpha[i] = np.pi/2

    Alpha = np.abs(Alpha) #np.abs(np.rad2deg(Alpha[W_fine[:-1]]))
    #print('alpha:',np.rad2deg(Alpha))
    Beta = np.abs(np.pi/2 - Alpha)

    dist = np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))
    area = np.pi*(r1+r2)*dist
    

    ##############################################
    #get W/s sputtered by D, C flux to wall
    ##############################################
     
    #get incoming ion energy and angle estimations where the integer input is z
    energyD, angleD = get_incoming_IEADs(1, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC1, angleC1 = get_incoming_IEADs(1, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC2, angleC2 = get_incoming_IEADs(2, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC3, angleC3 = get_incoming_IEADs(3, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC4, angleC4 = get_incoming_IEADs(4, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC5, angleC5 = get_incoming_IEADs(5, profiles, surfW, rmrsCoarse, rmrsFine)
    energyC6, angleC6 = get_incoming_IEADs(6, profiles, surfW, rmrsCoarse, rmrsFine)
    
    if plot_variables == 1:
        plt.close()
        plt.plot(rmrsFine, energyD, 'black', label='D1+')
        plt.plot(rmrsFine, energyC1, 'firebrick', label='C1+')
        plt.plot(rmrsFine, energyC2, 'darkorange', label='C2+')
        plt.plot(rmrsFine, energyC3, 'gold', label='C3+')
        plt.plot(rmrsFine, energyC4, 'limegreen', label='C4+')
        plt.plot(rmrsFine, energyC5, 'dodgerblue', label='C5+')
        plt.plot(rmrsFine, energyC6, 'mediumpurple', label='C6+')
        plt.plot(rmrsFine, 45.3362*np.ones(len(rmrsFine)), 'gray', label='W Eth')
        plt.legend()
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('energy [eV]')
        plt.title('Estimate of incident background ion energies')
        plt.savefig('plots/particle-source/incident_energy')
        
        plt.close()
        plt.plot(rmrsFine, angleD, 'black', label='D1+')
        plt.plot(rmrsFine, angleC1, 'firebrick', label='C1+')
        plt.plot(rmrsFine, angleC2, 'darkorange', label='C2+')
        plt.plot(rmrsFine, angleC3, 'gold', label='C3+')
        plt.plot(rmrsFine, angleC4, 'limegreen', label='C4+')
        plt.plot(rmrsFine, angleC5, 'dodgerblue', label='C5+')
        plt.plot(rmrsFine, angleC6, 'mediumpurple', label='C6+')        
        plt.legend()
        plt.xlabel('z [m]')
        plt.ylabel('angle [degrees]')
        plt.title('Used angles of incidence for incident background ions')
        plt.savefig('plots/particle-source/incident_angle.png')
        plt.close()

    #get sputtering yields for D0 and D1+ on W from fractal tridyn tables
    #Cspyld = get_ft_spyld(CsurfE, CsurfA, ftCFile)[0]
    spyldD = get_ft_spyld(1, energyD, angleD, ftDFile) #file input includes He, which we aren't using
    spyldC1 = get_ft_spyld(0, energyC1, angleC1, ftCFile)
    spyldC2 = get_ft_spyld(0, energyC2, angleC2, ftCFile)
    spyldC3 = get_ft_spyld(0, energyC3, angleC3, ftCFile)
    spyldC4 = get_ft_spyld(0, energyC4, angleC4, ftCFile)
    spyldC5 = get_ft_spyld(0, energyC5, angleC5, ftCFile)
    spyldC6 = get_ft_spyld(0, energyC6, angleC6, ftCFile)
    spyldW1 = get_ft_spyld(0, energyC1, angleC1, ftWFile)
    spyldW2 = get_ft_spyld(0, energyC2, angleC2, ftWFile)
    
    #get coarse flux profile from background D, C and refine to rmrsFine
    fluxCoarseD0 = np.abs(profiles.variables['flux_inner_target'][0][surfW])
    fluxCoarseD = np.abs(profiles.variables['flux_inner_target'][1][surfW])
    fluxCoarseC0 = np.abs(profiles.variables['flux_inner_target'][2][surfW])
    fluxCoarseC1 = np.abs(profiles.variables['flux_inner_target'][3][surfW])
    fluxCoarseC2 = np.abs(profiles.variables['flux_inner_target'][4][surfW])
    fluxCoarseC3 = np.abs(profiles.variables['flux_inner_target'][5][surfW])
    fluxCoarseC4 = np.abs(profiles.variables['flux_inner_target'][6][surfW])
    fluxCoarseC5 = np.abs(profiles.variables['flux_inner_target'][7][surfW])
    fluxCoarseC6 = np.abs(profiles.variables['flux_inner_target'][8][surfW])
    
    ffluxD0 = scii.interp1d(rmrsCoarse,fluxCoarseD0,fill_value='extrapolate')
    ffluxD = scii.interp1d(rmrsCoarse,fluxCoarseD,fill_value='extrapolate')
    ffluxC0 = scii.interp1d(rmrsCoarse,fluxCoarseC0,fill_value='extrapolate')
    ffluxC1 = scii.interp1d(rmrsCoarse,fluxCoarseC1,fill_value='extrapolate')
    ffluxC2 = scii.interp1d(rmrsCoarse,fluxCoarseC2,fill_value='extrapolate')
    ffluxC3 = scii.interp1d(rmrsCoarse,fluxCoarseC3,fill_value='extrapolate')
    ffluxC4 = scii.interp1d(rmrsCoarse,fluxCoarseC4,fill_value='extrapolate')
    ffluxC5 = scii.interp1d(rmrsCoarse,fluxCoarseC5,fill_value='extrapolate')
    ffluxC6 = scii.interp1d(rmrsCoarse,fluxCoarseC6,fill_value='extrapolate')
    
    fluxD0 = ffluxD0(rmrsFine)
    fluxD = ffluxD(rmrsFine)
    fluxC0 = ffluxC0(rmrsFine)
    fluxC1 = ffluxC1(rmrsFine)
    fluxC2 = ffluxC2(rmrsFine)
    fluxC3 = ffluxC3(rmrsFine)
    fluxC4 = ffluxC4(rmrsFine)
    fluxC5 = ffluxC5(rmrsFine)
    fluxC6 = ffluxC6(rmrsFine)
    fluxC = fluxC1 + fluxC2 + fluxC3 + fluxC4 + fluxC5 + fluxC6
    totalflux = fluxD0 + fluxD + fluxC0 + fluxC
    Cfraction = np.sum(fluxC0 + fluxC) / np.sum(totalflux)
    print('C Fraction of Total Incoming Flux:', Cfraction)

    #multiply incoming ion flux by Y_s to get sputtered W flux by each species
    sputt_fluxD = spyldD*fluxD
    sputt_fluxC1 = spyldC1*fluxC1
    sputt_fluxC2 = spyldC2*fluxC2
    sputt_fluxC3 = spyldC3*fluxC3
    sputt_fluxC4 = spyldC4*fluxC4
    sputt_fluxC5 = spyldC5*fluxC5
    sputt_fluxC6 = spyldC6*fluxC6
    sputt_flux = sputt_fluxD + sputt_fluxC1 + sputt_fluxC2 + sputt_fluxC3 + sputt_fluxC4 + sputt_fluxC5 + sputt_fluxC6
    print('SPUTT FLUX',len(sputt_flux),'\n',sputt_flux)
    print('\n')
    print('W eroded flux per nP:', np.sum(sputt_flux)/nP, 'm-2 s-1')


    #multiply by area to get the outgoing particles per second
    pps = np.multiply(sputt_flux,area)
    print('Total actual W eroded per second:', np.sum(pps), 's-1 \n')
    pps_weights = nP*pps/np.sum(pps)

    if plot_variables == 1: 
        plt.rcParams.update({'font.size':14})
        plt.close()
        if tile_shift_indices != []:
            for i,v in enumerate(tile_shift_indices):
                if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
                else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
        if Bangle_shift_indices != []:
            for i,v in enumerate(Bangle_shift_indices):
                if i==0: plt.axvline(x=rmrsCoarse[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
                else: plt.axvline(x=rmrsCoarse[v], color='k', linestyle='dotted')
        
        plt.plot(rmrsFine, fluxC1, 'firebrick', label='C$^{1+}$')
        plt.plot(rmrsFine, fluxC2, 'darkorange', label='C$^{2+}$')
        plt.plot(rmrsFine, fluxC3, 'gold', label='C$^{3+}$')
        plt.plot(rmrsFine, fluxC4, 'limegreen', label='C$^{4+}$')
        plt.plot(rmrsFine, fluxC5, 'dodgerblue', label='C$^{5+}$')
        plt.plot(rmrsFine, fluxC6, 'mediumpurple', label='C$^{6+}$')
        '''
        plt.plot(rmrsCoarse, fluxCoarseC1, 'firebrick', linestyle='dotted', label='C$^{1+}$')
        plt.plot(rmrsCoarse, fluxCoarseC2, 'darkorange', linestyle='dotted', label='C$^{2+}$')
        plt.plot(rmrsCoarse, fluxCoarseC3, 'gold', linestyle='dotted', label='C$^{3+}$')
        plt.plot(rmrsCoarse, fluxCoarseC4, 'limegreen', linestyle='dotted', label='C$^{4+}$')
        plt.plot(rmrsCoarse, fluxCoarseC5, 'dodgerblue', linestyle='dotted', label='C$^{5+}$')
        plt.plot(rmrsCoarse, fluxCoarseC6, 'mediumpurple', linestyle='dotted', label='C$^{6+}$')
        '''
        plt.xlim(right=0.14)
        plt.ylim(top=2.05e20)
        plt.yticks(1e20*np.array([0,0.5, 1, 1.5, 2.0]))
        #plt.yscale('log')
        plt.xlabel('D-Dsep [m]', fontsize=14)
        plt.ylabel('Flux [m$^{-2}$s$^{-1}$]', fontsize=16)
        plt.legend(loc='upper right', fontsize=12)
        plt.title('C flux to W surface', fontsize=24)
        plt.show(block=False)
        plt.savefig('plots/particle-source/incident_flux.png')

        plt.close()
        #plt.rcParams.update({'font.size':16})
        #plt.rcParams.update({'lines.linewidth':3})
        
        plt.plot(rmrsFine, spyldD, 'black', label='D$^{1+}$')
        plt.plot(rmrsFine, spyldC1, 'firebrick', label='C$^{1+}$')
        plt.plot(rmrsFine, spyldC2, 'darkorange', label='C$^{2+}$')
        plt.plot(rmrsFine, spyldC3, 'gold', label='C$^{3+}$')
        plt.plot(rmrsFine, spyldC4, 'limegreen', label='C$^{4+}$')
        plt.plot(rmrsFine, spyldC5, 'dodgerblue', label='C$^{5+}$')
        plt.plot(rmrsFine, spyldC6, 'mediumpurple', label='C$^{6+}$')
        plt.plot(rmrsFine, spyldW1, 'rosybrown', label='W$^{1+}$')
        plt.plot(rmrsFine, spyldW2, 'burlywood', label='W$^{2+}$')
        #plt.xlim(-0.025)
        plt.legend(fontsize=12)
        plt.xlabel('D-Dsep [m]', fontsize=14)
        plt.ylabel('Sputtering Yield', fontsize=16)
        plt.title('W sputtering yield by incident ions',fontsize=24)
        plt.show(block=False)
        plt.savefig('plots/particle-source/spyld.png')
        
        plt.close()
        if tile_shift_indices != []:
            for i,v in enumerate(tile_shift_indices):
                if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Wall\nVertices')
                else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
        if Bangle_shift_indices != []:
            for i,v in enumerate(Bangle_shift_indices):
                if i==0: plt.axvline(x=rmrsCoarse[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
                else: plt.axvline(x=rmrsCoarse[v], color='k', linestyle='dotted')
        
        plt.plot(rmrsFine, sputt_fluxD, 'black', label='D$^{1+}$')
        plt.plot(rmrsFine, sputt_fluxC1, 'firebrick', label='C$^{1+}$')
        plt.plot(rmrsFine, sputt_fluxC2, 'darkorange', label='C$^{2+}$')
        plt.plot(rmrsFine, sputt_fluxC3, 'gold', label='C$^{3+}$')
        plt.plot(rmrsFine, sputt_fluxC4, 'limegreen', label='C$^{4+}$')
        plt.plot(rmrsFine, sputt_fluxC5, 'dodgerblue', label='C$^{5+}$')
        plt.plot(rmrsFine, sputt_fluxC6, 'mediumpurple', label='C$^{6+}$')
        #plt.xlim(-0.05)
        plt.xlabel('D-Dsep [m]', fontsize=14)
        plt.ylabel('Flux [m$^{-2}$s$^{-1}$]', fontsize=16)
        plt.legend(loc='upper left', fontsize=12)
        plt.title('Flux of sputtered W', fontsize=24)
        plt.show(block=False)
        plt.savefig('plots/particle-source/sputt_flux_charge_dependent.png')

        plt.close()
        plt.plot(rmrsFine, sputt_flux)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Flux [#/m2s]')
        plt.title('Flux of W sputtered off wall')
        plt.savefig('plots/particle-source/sputt_flux.png') 

        plt.close()
        area_cm = area*100*100
        plt.plot(rmrsFine, area_cm)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Area [cm2]')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(1,5))        
        plt.title('Toroidal surface area of each line segment')
        plt.savefig('plots/particle-source/area.png')
        
        plt.close()
        plt.plot(rmrsFine, pps_weights)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Weighted PPS [#/s]')
        plt.title('Computational Weights for Initially Sputtered W/s \n nP = '+str(nP))
        plt.savefig('plots/particle-source/pps_weights.png')

    #########################################
    #get x,y,z distributions for sputtered W
    #########################################

    #confirm nP stays constant
    for i in range(len(pps_weights)): pps_weights[i] = round(pps_weights[i])
    nP_diff = int(nP-np.sum(pps_weights))
    
    if nP_diff > 0:
        for i in range(abs(nP_diff)):
            rand_index = np.random.choice(len(pps_weights))
            pps_weights[rand_index] += 1
    elif nP_diff < 0:
        while nP_diff<0:
            rand_index = np.random.choice(len(pps_weights))
            if pps_weights[rand_index] - 1 > 0:
                pps_weights[rand_index] -= 1
                nP_diff+=1
    
    pps_weights = np.round(pps_weights)
    int_weights = np.array(pps_weights,dtype='int')
    nP_diff = nP-np.sum(int_weights)
    pps_per_nP = np.sum(pps)/nP

    print('total nP', nP)
    print('pps over nP', pps_per_nP)
    print('nP(r_mid):', int_weights)
    print('nP_diff should be 0: ', nP_diff)

    #define adjustment into the sheath because particles can't start exactly on the wall
    adj = 2e-7

    #populate x,y,z with r_mid,0,z_mid
    if configuration == 'random': 
        x,y,z = random(nP,int_weights,adj,slope,Beta, r1,z1,r2,z2)
    elif configuration == 'uniform': 
        x,y,z = uniform(nP,int_weights,adj,slope,Beta, r1,z1,r2,z2)
    elif configuration == 'midpoint': 
        x,y,z = midpoints(nP,int_weights, adj,slope,Beta, r1,z1,r2,z2)
    else:
        print('(x,y,z) configuration not set')

    #########################################
    #get vx,vy,vz from IEADs
    #########################################

    #use PyGITR to set up vx,vy,vz,E,theta,psi distributions
    PartDist = Particles.ParticleDistribution(nP, ListAttr=['vx','vy','vz'])

    vx = np.empty(0)
    vy = np.empty(0)
    vz = np.empty(0)
    
    z_coarse_index = 0 #debugging
    #print('COARSE SEGMENT: 0') #debugging
    
    for i in range(len(int_weights)): #range(29,34):#
    
        weight = int(int_weights[i])
        weight_diff = np.zeros(len(int_weights))
        m = np.sign(slope[i])

        if weight>0:
            
            if use_fractal_tridyn_outgoing_IEADS:
                
                E = np.empty(0)
                PolAng = np.empty(0)
                AziAng = np.empty(0)
                
                #confirm species-weighted pps integer weights stay constant
                weightsC1 = weight * sputt_fluxC1[i]/sputt_flux[i]
                weightsC2 = weight * sputt_fluxC2[i]/sputt_flux[i]
                weightsC3 = weight * sputt_fluxC3[i]/sputt_flux[i]
                weightsC4 = weight * sputt_fluxC4[i]/sputt_flux[i]
                weightsC5 = weight * sputt_fluxC5[i]/sputt_flux[i]
                weightsC6 = weight * sputt_fluxC6[i]/sputt_flux[i]
                species_weights = np.array([weightsC1,weightsC2,weightsC3,weightsC4,weightsC5,weightsC6])
                
                for j in range(len(species_weights)): 
                    species_weights[j] = int(round(species_weights[j]))
                
                #randomly add particles to species weight to adjust for rounding                    
                weight_diff[i] = round(weight - np.sum(species_weights))
                
                if weight_diff[i] > 0:
                    for d in range(abs(int(weight_diff[i]))):
                        rand_index = np.random.choice(len(species_weights))
                        species_weights[rand_index] += 1
                elif weight_diff[i] < 0:
                    for d in range(abs(int(weight_diff[i]))):
                        pos_indices = np.where(species_weights>0)
                        rand_index = np.random.choice(pos_indices[0])
                        species_weights[rand_index] -= 1
                
                for j in range(len(species_weights)): 
                    species_weights[j] = round(species_weights[j])
                weight_diff[i] = int(weight - np.sum(species_weights))

                weightC1 = int(species_weights[0])
                weightC2 = int(species_weights[1])
                weightC3 = int(species_weights[2])
                weightC4 = int(species_weights[3])
                weightC5 = int(species_weights[4])
                weightC6 = int(species_weights[5])
                
                #get IEADs for outgoing W from from fractal tridyn
                EC1, PhiC1, ThetaC1 = get_ft_IEAD(energyC1[i], weightC1, ftCFile)
                EC2, PhiC2, ThetaC2 = get_ft_IEAD(energyC2[i], weightC2, ftCFile)
                EC3, PhiC3, ThetaC3 = get_ft_IEAD(energyC3[i], weightC3, ftCFile)
                EC4, PhiC4, ThetaC4 = get_ft_IEAD(energyC4[i], weightC4, ftCFile)
                EC5, PhiC5, ThetaC5 = get_ft_IEAD(energyC5[i], weightC5, ftCFile)
                EC6, PhiC6, ThetaC6 = get_ft_IEAD(energyC6[i], weightC6, ftCFile)
                
                #populate 1D list of energies, polar angles, and azimuthal angles
                for w in range(weightC1):
                    E = np.append(E,EC1[w])
                    PolAng = np.append(PolAng,PhiC1[w])
                    AziAng = np.append(AziAng,ThetaC1[w])
                for w in range(weightC2):
                    E = np.append(E,EC2[w])
                    PolAng = np.append(PolAng,PhiC2[w])
                    AziAng = np.append(AziAng,ThetaC2[w])
                for w in range(weightC3):
                    E = np.append(E,EC3[w])
                    PolAng = np.append(PolAng,PhiC3[w])
                    AziAng = np.append(AziAng,ThetaC3[w])
                for w in range(weightC4):
                    E = np.append(E,EC4[w])
                    PolAng = np.append(PolAng,PhiC4[w])
                    AziAng = np.append(AziAng,ThetaC4[w])
                for w in range(weightC5):
                    E = np.append(E,EC5[w])
                    PolAng = np.append(PolAng,PhiC5[w])
                    AziAng = np.append(AziAng,ThetaC5[w])
                for w in range(weightC6):
                    E = np.append(E,EC6[w])
                    PolAng = np.append(PolAng,PhiC6[w])
                    AziAng = np.append(AziAng,ThetaC6[w])
                
                PolAng = np.deg2rad(PolAng)
                AziAng = np.deg2rad(AziAng)
                
                #debugging
                '''
                #convert IEADs to vx,vy,vz unit vectors in particle frame of ref
                vx_prime = np.multiply(np.cos(PolAng), np.cos(AziAng))
                vy_prime = np.multiply(np.cos(PolAng), np.sin(AziAng))
                vz_prime = m*np.sin(PolAng)
                (vx_lab,vy_lab,vz_lab) = (vx_prime,vy_prime,vz_prime)
                
                z_inner_target = z_right_target[surfW]
                if z2[i] < z_inner_target[z_coarse_index+1]:
                    z_coarse_index += 1
                    print('\n COARSE SEGMENT:', z_coarse_index, i)
                
                
                #make histograms of energies for each line segment
                plt.close()
                plt.hist(E,50) #in eV
                plt.xlabel('Energy [eV]')
                plt.ylabel('Counts')
                plt.title('Histogram of Initial Outgoing Energies \n Segment: '+str(i))
                plt.show(block=True)
                
                #make histograms of polar angles for each line segment
                plt.close()
                plt.hist(PolAng, 50) #in degrees
                plt.xlabel('Polar Angle [degree]')
                plt.ylabel('Counts')
                plt.title('Histogram of Initial Outgoing Polar Angles \n Segment: '+str(i))
                plt.show(block=True)
                '''
        
            else:
            
                #get IEADs for sputtered W
                E = PartDist.Generate(weight, 'Thompson')
                PolAng = PartDist.Generate(weight, 'SinCos', x=np.linspace(0,np.pi/2,10*weight))
                AziAng = PartDist.Generate(weight, 'Gaussian', x=np.linspace(-np.pi,np.pi,10*weight))
            
            #convert IEADs to vx,vy,vz unit vectors in particle frame of ref
            vx_prime = -np.cos(PolAng)
            vy_prime = np.multiply(np.sin(PolAng), -np.cos(AziAng))
            vz_prime = np.multiply(np.sin(PolAng), -np.sin(AziAng))
            PartDist.SetAttr('vx', vx_prime)
            PartDist.SetAttr('vy', vy_prime)
            PartDist.SetAttr('vz', vz_prime)
    
            #rotate vx,vy,vz from particle frame to lab frame
            PartDist.RotateAngle('v', m*Beta[i],0, Degree=False)
            
            vx_lab = PartDist.Particles['vx']
            vy_lab = vy_prime
            vz_lab = PartDist.Particles['vz']
    
            #convert unit vectors to vx,vy,vz
            W_kg = 183.84 * 1.6605e-27 #mass of W in kg
            vtot = np.sqrt(2*E*1.6022e-19/W_kg) #convert eV to m/s
            
            vx = np.append(vx, vtot*vx_lab)
            vy = np.append(vy, vtot*vy_lab)
            vz = np.append(vz, vtot*vz_lab)
            
            '''#debugging
            #print(vz_lab)
            plt.close()
            plt.scatter(vx_lab,vz_lab,s=1)
            plt.axis('Scaled')
            plt.xlabel('vx')
            plt.ylabel('vz')
            plt.title('Polar Angle Distribution in the Lab Frame \n nP='+str(int_weights[i]))
            plt.show(block=True)'''
    '''
    print('vx',np.average(vx))
    print('vy',np.average(vy))
    print('vz',np.average(vz))
    
    #debugging
    blocker=False
    print(np.cumsum(int_weights))
    first_particle = np.cumsum(int_weights)[28] + 1
    last_particle = np.cumsum(int_weights)[33] - first_particle + 1
    print('first particle:', first_particle)
    print('last particle:', last_particle)
    
    print('x:',x[first_particle:last_particle])
    print('y:',y[first_particle:last_particle])
    print('z:',z[first_particle:last_particle])
    plt.close()
    plt.hist(vx[first_particle:last_particle])
    plt.title('vx')
    plt.show(block=blocker)
    plt.hist(vy[first_particle:last_particle])
    plt.title('vy')
    plt.show(block=blocker)
    plt.hist(vz[first_particle:last_particle])
    plt.title('vz')
    plt.show(block=blocker)'''

    #double check that all particles actually received an energy and 2 angles
    print('species weights_diff should be 0:', np.sum(weight_diff))
    
    if plot_variables == 1:
        #plot Thomson E dist
        plt.close()
        plt.hist(E,bins=54)
        plt.xlabel('Energy Bins [eV]')
        plt.ylabel('Histogram')
        plt.title('Thomson Energy Distribution \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/thomsonDist.png')
        
        #plot particle framed v_dist relations
        plt.close()
        plt.scatter(vx_prime,vy_prime,s=1)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vy')
        plt.title('Azimuthal Angle in the Particle Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vxvy_prime.png')
        
        plt.close()
        plt.scatter(vx_prime,vz_prime,s=1)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vz')
        plt.title('Azimuthal Angle in the Particle Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vxvz_prime.png')
        
        plt.close()
        plt.scatter(vy_prime,vz_prime,s=1)
        plt.axis('Scaled')
        plt.xlabel('vy')
        plt.ylabel('vz')
        plt.title('Polar Angle in the Particle Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vyvz_prime.png')
        
        plt.close()
        plt.scatter(vx_lab,vy_lab,s=1)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vy')
        plt.title('Azimuthal Angle in the Lab Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vxvy_lab.png')
        
        plt.close()
        plt.scatter(vx_lab,vz_lab,s=1)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vz')
        plt.title('Polar Angle Distribution in the Lab Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vxvz_lab.png')
        
        plt.close()
        plt.scatter(vy_lab,vz_lab,s=1)
        plt.axis('Scaled')
        plt.xlabel('vy')
        plt.ylabel('vz')
        plt.title('Azimuthal Angle in the Lab Frame \n nP='+str(int_weights[-1]))
        plt.savefig('plots/particle-source/vyvz_lab.png')
        plt.close()
        
        

    #########################################
    #make NetCDF Particle Source file
    #########################################
    '''
    #debugging
    nP = last_particle - first_particle
    print('debugging nP:',nP)
    x = x[first_particle:last_particle]
    y = y[first_particle:last_particle]
    z = z[first_particle:last_particle]
    vx = vx[first_particle:last_particle]
    vy = vy[first_particle:last_particle]
    vz = vz[first_particle:last_particle]
    '''
    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nP)
    onee = rootgrp.createDimension("one",1)
    PPS_per_nP = rootgrp.createVariable("pps_per_nP","f8",("one"))
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    PPS_per_nP[:] = pps_per_nP
    xxx[:] = x
    yyy[:] = y
    zzz[:] = z
    vxx[:] = vx
    vyy[:] = vy
    vzz[:] = vz
    rootgrp.close()
    
    print('Created particleSource.nc')
    
    return pps_per_nP, sputt_flux, fluxD, fluxC



def midpoints(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
    #get midpoints of coords
    r_mid = np.zeros(len(r1)-1)
    z_mid = np.zeros(len(z1)-1)
    for i in range(len(r1)):
        r_mid[i] = np.average(np.array([r1[i],r2[i]]))
        z_mid[i] = np.average(np.array([z1[i],z2[i]]))
        
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        x[counter:counter+pps_weights[i]] = r_mid[i] - adj*np.abs(np.cos(Beta[i]))
        z[counter:counter+pps_weights[i]] = z_mid[i] + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]
        
    return x,y,z

def uniform(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    tally = 0
    for i in range(len(pps_weights)):
        if pps_weights[i]>0:
            dr = (r2[i]-r1[i])/pps_weights[i]
            dz = (z2[i]-z1[i])/pps_weights[i]
            r_segment = np.linspace(r1[i]+dr/2, r2[i]-dr/2, pps_weights[i])
            z_segment = np.linspace(z1[i]+dz/2, z2[i]-dz/2, pps_weights[i])
            
            x[tally:tally+pps_weights[i]] = r_segment - adj*np.abs(np.cos(Beta[i]))
            z[tally:tally+pps_weights[i]] = z_segment + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
            tally += pps_weights[i]
            
    return x,y,z

def random(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        for j in range(pps_weights[i]):
            chi = np.random.rand(1)
            x[counter+j] = r1[i]+chi*(r2[i]-r1[i]) - adj*np.abs(np.cos(Beta[i]))
            z[counter+j] = z1[i]+chi*(z2[i]-z1[i]) + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]
        
    return x,y,z



def refine_linear(a,b,coarse_rmrs, fine_rmrs, coarse_variable):
    fine_variable = np.zeros(len(fine_rmrs))
    for i in range(len(fine_rmrs)):
        #find the smallest coarse index near the fine index
        index = np.where(coarse_rmrs<=fine_rmrs[i])
        if index[0].any() > 0:
            lower_bound_index = np.max(index)
        else: lower_bound_index = 0
        
        if lower_bound_index < len(coarse_rmrs)-1:
            upper_bound_index = lower_bound_index + 1 
            lower_bound_rmrs = coarse_rmrs[lower_bound_index]
            upper_bound_rmrs = coarse_rmrs[upper_bound_index]
            lower_bound_variable = coarse_variable[lower_bound_index]
            upper_bound_variable = coarse_variable[upper_bound_index]
            
            fraction = (fine_rmrs[i]-lower_bound_rmrs) / (upper_bound_rmrs-lower_bound_rmrs)
            fine_variable[i] = fraction * (upper_bound_variable-lower_bound_variable) + lower_bound_variable
        else: fine_variable[i] = coarse_variable[-1]
        
    return fine_variable

def get_incoming_IEADs(q, profiles, surfW, rmrsCoarse, rmrsFine):
    #extract plasma parameters at the wall indices
    teCoarse = profiles.variables['te_inner_target'][surfW]
    tiCoarse = profiles.variables['ti_inner_target'][surfW]
    
    #get temp as a function of rmrsCoarse
    fte = scii.interp1d(rmrsCoarse, teCoarse, fill_value='extrapolate')
    fti = scii.interp1d(rmrsCoarse, tiCoarse, fill_value='extrapolate')
    
    #interpolate temp at points in rmrsFine
    te = fte(rmrsFine)
    ti = fti(rmrsFine)
    
    SimpleEnergyEst = 2*ti+3*te*q
    
    if q<0:
        Esp, f, b, c, ThetaMax = FittingParameters_NonW(SimpleEnergyEst)
        AngleEst = ThetaMax
    else: AngleEst = 70*np.ones(len(SimpleEnergyEst))
    
    return SimpleEnergyEst, AngleEst

def get_ft_spyld(S, surfE, surfA, ftBFile):
    #import sputtering yield tables for incident ions on W
    ftB = netCDF4.Dataset(ftBFile, "r", format="NETCDF4")
    ftE = ftB.variables['E'][:]
    ftA = ftB.variables['A'][:]
    
    #handle surface energies less than the minimum f-TRIDYN grid energy
    surfE[np.where(surfE<min(ftE))] = min(ftE) 
    
    try: 
        spyld = ftB.variables['spyld'][S][:]
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA))
    
    except: 
        spyld = ftB.variables['spyld'][:]
        surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA))
    
    return surfY

def get_ft_IEAD(surfE, weight, ftBFile):
    #CANNOT TAKE FULL SURFACES, SURFE AND SURFA MUST BE SINGLE VARIABLES
    #import IEAD tables for incident ions on W
    ftB = netCDF4.Dataset(ftBFile, "r", format="NETCDF4")
    
    #find closest ftE value to plug into distribution tables
    ftE = ftB.variables['E'][:]
    Ediff = 1000
    for i,v in enumerate(ftE):
        Ediff_test = np.abs(v-surfE)
        if Ediff_test<Ediff:
            Ediff=Ediff_test
            nearestEindex = i
                
    #find closest ftA value to plug into distribution tables
    nearestAindex = 4

    #get PDFs and CDFs over energy/angle grids
    energyGrid = ftB.variables['eDistEgrid'][:]
    energyPDF = ftB.variables['energyDist'][0][nearestEindex][nearestAindex][:]
    energyCDF = np.cumsum(energyPDF)/np.sum(energyPDF)
    
    phiGrid = ftB.variables['phiGrid'][:]
    cosXPDF = ftB.variables['cosXDist'][0][nearestEindex][nearestAindex][:]
    cosXCDF = np.cumsum(cosXPDF)/np.sum(cosXPDF)
    
    thetaGrid = ftB.variables['thetaGrid'][:]
    cosYPDF = ftB.variables['cosYDist'][0][nearestEindex][nearestAindex][:]
    cosYCDF = np.cumsum(cosYPDF)/np.sum(cosYPDF)
    
    #stretch thetaGrid to span 0 to 360
    thetaGrid = 2*thetaGrid-180+45
    
    if np.max(energyPDF)==0:
        energySampled = np.zeros(weight)
        phiSampled = np.zeros(weight)
        thetaSampled = np.zeros(weight)
    
    else:
        #turn cdfs into functions from which to randomly sample energies and angles 
        energyFunc = scii.interp1d(energyCDF, energyGrid, bounds_error=False, fill_value=(energyGrid[0], energyGrid[-1]))
        phiFunc = scii.interp1d(cosXCDF, phiGrid,  bounds_error=False, fill_value=(0, 90))
        thetaFunc = scii.interp1d(cosYCDF, thetaGrid, bounds_error=False, fill_value=(0, 180))
        
        #debugging
        '''
        plt.close()
        plt.plot(phiGrid,cosXPDF)
        plt.title('phi')
        plt.xlabel('degree')
        plt.ylabel('counts')
        plt.show(block=True)
        plt.plot(thetaGrid,cosYPDF)
        plt.title('theta')
        plt.xlabel('degree')
        plt.ylabel('counts')
        plt.show(block=True)
        '''
        
        #randomly sample energies and angles from CDFs
        [energyXi, phiXi, thetaXi] = np.random.rand(3, weight)
        energySampled = energyFunc(energyXi)
        phiSampled = phiFunc(phiXi)
        thetaSampled = thetaFunc(thetaXi)    
    
    return energySampled, phiSampled, thetaSampled

def FittingParameters_NonW(E):
    Eth=45.3362
    Esp = np.ones(len(E))
    f = np.zeros(len(E))
    b = np.zeros(len(E))
    c = np.zeros(len(E))
    ThetaMax = np.zeros(len(E))
    for i,v in enumerate(E):
        if v<Eth:
            continue
        elif v>=Eth and v<49:
            f[i],b[i],c[i] = 2.9557, 5.8879, 0.9465
        elif v>=49 and v<51:
            f[i],b[i],c[i] = 1.7735, 4.3144, 0.9468
        elif v>=51 and v<53.5:
            f[i],b[i],c[i] = 1.2707, 3.6458, 0.8840
        elif v>=53.5 and v<57.5:
            f[i],b[i],c[i] = 1.1002, 3.3751, 0.9604
        elif v>=57.5 and v<65:
            f[i],b[i],c[i] = 0.4622, 2.5095, 1.0118
        elif v>=65 and v<75:
            f[i],b[i],c[i] = 0.4622, 2.5095, 1.0118
        elif v>=75 and v<85:
            f[i],b[i],c[i] = 2.7960, 3.4029, 0.8841
        elif v>=85 and v<95:
            f[i],b[i],c[i] = 2.1152, 2.6541, 0.9226
        elif v>=95 and v<110:
            f[i],b[i],c[i] = 1.7312, 2.1735, 0.9489
        elif v>=110 and v<130:
            f[i],b[i],c[i] = 1.6230, 1.6737, 1.0004
        elif v>=130 and v<170:
            f[i],b[i],c[i] = 1.7195, 1.5092, 1.0176
            ThetaMax[i] = 30.54
        elif v>=170 and v<250:
            f[i],b[i],c[i] = 2.0138, 1.3460, 1.0316
            ThetaMax[i] = 51.98
        elif v>=250 and v<400:
            f[i],b[i],c[i] = 2.2531, 1.2151, 1.0310
            ThetaMax[i] = 59.47
        elif v>=400 and v<750:
            f[i],b[i],c[i] = 2.4324, 1.1313, 1.0171
            ThetaMax[i] = 62.20
        elif v>=750 and v<2000:
            f[i],b[i],c[i] = 2.4383, 0.9940, 0.9936
            ThetaMax[i] = 66.00
        else:
            print('WARNING: TABULAR SPUTTERING FITTING PARAMETERS MISSING')
    return Esp, f, b, c, ThetaMax

def get_analytic_spyld(surfE, surfA, Z1=6, M1=12, Z2=74, M2=183.84, \
                       FitParam='N', Eth=45.3362, lam=0.0921, q=1.4389, mu=2.0225):
    #this entire function that is simply the Eckstein formula
    #defaults are for C on W but with fitting parameters, Eth, and Esp for N on W
    #energies are in eV and angles are in degrees
    #M1, Z1 are for the projectile
    #M2, Z2 are for the target
    #Eth is the threshold energy for any sputtering to occur
    #Esp is the SBE for self-bombardment
    #lam, q, mu, f, b, c are fitting parameters from Eckstein tables
    
    e = 14.399651
    
    a_L = 0.0529177*((9*np.pi**2/128)**(1/3)) + \
        (Z1**(2/3) + Z2**(2/3))**-0.5
        
    epsilon_L = surfE * (M2/(M1+M2)) * (a_L/(Z1*Z2*e**2))
    
    omega = epsilon_L + 0.1728*np.sqrt(epsilon_L) + \
        0.008*epsilon_L**0.1504
        
    snKrC = 0.5*np.log(1+1.2288*epsilon_L)/omega
    
    #spyld(E, angle = 0) normal incidence
    Y_0 = np.zeros(len(surfE))
    for i,v in enumerate(surfE):
        if v>=Eth:
            Y_0[i] = q*snKrC[i]*(v/Eth-1)**mu / \
                (lam/omega[i] + (v/Eth-1)**mu)
    
    #choose basis for fitting parameters for different energies
    Esp, f, b, c, ThetaMax = FittingParameters_NonW(surfE)
    
    theta_0 = np.pi - np.arccos(np.sqrt(1/(1+surfE/Esp)))
    
    #spyld(E,A)    
    Y = Y_0 * (np.cos((np.pi*surfA/(2*theta_0))**c))**(-f) * \
        np.exp(b*(1-1/np.cos((np.pi*surfA/(2*theta_0))**c)))
        
    return Y_0



if __name__ == "__main__":
    #init()
    
    distributed_source(nP=int(2e3), surfW=np.arange(16,25), \
                tile_shift_indices = [2,6], \
                Bangle_shift_indices = [3,6], \
                geom = '../input/gitrGeometry.cfg', \
                profiles_file = '../input/plasmaProfiles.nc', \
                gitr_rz = 'assets/gitr_rz.txt', \
                rmrs_fine_file = 'assets/rmrs_fine.txt', \
                W_fine_file = 'assets/W_fine.txt', \
                ftDFile = 'assets/ftridynBackgroundD.nc', \
                ftCFile = 'assets/ftridynBackgroundC.nc', \
                configuration = 'random', \
                plot_variables = 1)



