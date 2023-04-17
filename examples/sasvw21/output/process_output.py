import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import solps

def init():
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    
    #check which target is the region of interest
    W_indices = np.arange(16,31)
    print(W_indices)
    r_inner_target = profiles.variables['r_inner_target'][:]
    z_inner_target = profiles.variables['z_inner_target'][:]
    rmrs = profiles.variables['rmrs_inner_target_midpoints'][W_indices]
    
    #plot r,z to check that target indices are restricted to W surfaces
    plt.close()
    plt.plot(r_inner_target, z_inner_target)
    
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':1.2})
    plt.rcParams.update({'lines.markersize':1})

    return profiles, W_indices, r_inner_target, z_inner_target, rmrs
profiles, W_indices, R, Z, rmrs = init()

def plot_gitr_gridspace():
    #import vars from plasmaProfiles.nc
    gridr = profiles.variables['gridr'][:]
    gridz = profiles.variables['gridz'][:]
    values = profiles.variables['values'][:]
    te = profiles.variables['te'][:]
    
    #set off grid space to 0
    off_grid_inds = np.where(te <= 0.0)
    te[off_grid_inds] = np.nan;
    
    #reformat
    nx,ny,crx,cry,region = solps.read_b2f_geometry('../setup/assets/b2fgmtry')
    meshgridr,meshgridz = np.meshgrid(gridr,gridz)
    print(meshgridr.shape)
    
    #plot gridspace
    plt.close()
    plt.plot(R,Z, '-k', linewidth=2)
    plt.pcolormesh(gridr,gridz,te,shading='nearest')
    plt.colorbar()
    plt.xlim(1.4,1.55)
    plt.ylim(1.05,1.25)
    return
    #store plasma parameters for cells along the W surfaces as a function of rmrs
    Bmag = profiles.variables['Bmag_inner_target'][W_indices]
    Bangle = profiles.variables['Bangle_inner_target'][W_indices]
    te = profiles.variables['te_inner_target'][W_indices]
    ti = profiles.variables['ti_inner_target'][W_indices]
    ne = profiles.variables['ne_inner_target'][W_indices]
    ni = profiles.variables['ni_inner_target'][:,W_indices]
    flux = profiles.variables['flux_inner_target'][:,W_indices]
    
    #plot magnitude of B-field
    plt.close()
    plt.plot(rmrs,Bmag)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Magnetic Flux Density [T]')
    plt.title('B-field magnitude along W surface')
    plt.savefig('plots/surface_profiles/Bmag.png')    
    
    #plot angle between B-field and normal incidence to the wall
    plt.close()
    plt.plot(rmrs,Bangle)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Angle [$\circ$]')
    plt.title('B-field angle along W surface')
    plt.savefig('plots/surface_profiles/Bangle.png')
    
    #plot te along the W surfaces
    plt.close()
    plt.plot(rmrs,te)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Temperature [eV]')
    plt.title('Electron temperature along W surface')
    plt.savefig('plots/surface_profiles/te.png')
    
    #plot ti along the W surfaces
    plt.close()
    plt.plot(rmrs,ti)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Temperature [eV]')
    plt.title('Ion temperature along W surface')
    plt.savefig('plots/surface_profiles/ti.png')
    
    #plot ne along the W surfaces
    plt.close()
    plt.plot(rmrs,ne)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Density [m$^{-3}$]')
    plt.title('Electron density along W surface')
    plt.savefig('plots/surface_profiles/ne.png')
    
    #plot average ni along the W surfaces
    plt.close()
    ni_avg = np.average(ni, axis=0)
    plt.plot(rmrs,ni_avg)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Density [m$^{-3}$]')
    plt.title('Average ion density along W surface')
    plt.savefig('plots/surface_profiles/niAvg.png')
    
    #plot ni along the W surfaces for deuterium
    plt.close()
    plt.plot(rmrs,ni[0],color='k',label='D 0')
    plt.plot(rmrs,ni[1],color='firebrick',label='D 1+')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Density [m$^{-3}$]')
    plt.title('D ion density along W surface')
    plt.legend()
    plt.savefig('plots/surface_profiles/niD.png')
    
    #plot flux along the W surfaces for deuterium
    plt.close()
    plt.plot(rmrs,flux[0],color='k',label='D 0')
    plt.plot(rmrs,flux[1],color='firebrick',label='D 1+')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [m$^{-2}$s$^{-1}$]')
    plt.title('D ion flux along W surface')
    plt.legend()
    plt.savefig('plots/surface_profiles/fluxD.png')
    
    #plot ni along the W surfaces for carbon
    plt.close()
    plt.plot(rmrs,ni[2],color='k',label='C 0')
    plt.plot(rmrs,ni[3],color='firebrick',label='C 1+')
    plt.plot(rmrs,ni[4],color='darkorange',label='C 2+')
    plt.plot(rmrs,ni[5],color='gold',label='C 3+')
    plt.plot(rmrs,ni[6],color='limegreen',label='C 4+')
    plt.plot(rmrs,ni[7],color='dodgerblue',label='C 5+')
    plt.plot(rmrs,ni[8],color='mediumslateblue',label='C 6+')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Density [m$^{-3}$]')
    plt.title('C ion density along W surface')
    plt.legend(loc='lower right')
    plt.savefig('plots/surface_profiles/niC.png')
    
    #plot flux along the W surfaces for carbon
    plt.close()
    plt.plot(rmrs,flux[2],color='k',label='C 0')
    plt.plot(rmrs,flux[3],color='firebrick',label='C 1+')
    plt.plot(rmrs,flux[4],color='darkorange',label='C 2+')
    plt.plot(rmrs,flux[5],color='gold',label='C 3+')
    plt.plot(rmrs,flux[6],color='limegreen',label='C 4+')
    plt.plot(rmrs,ni[7],color='dodgerblue',label='C 5+')
    plt.plot(rmrs,ni[8],color='mediumslateblue',label='C 6+')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [m$^{-2}$s$^{-1}$]')
    plt.yticks(fontsize=10.5)
    plt.title('C ion flux along W surface')
    plt.legend()
    plt.savefig('plots/surface_profiles/fluxC.png')
    
    #electron temperature dependence of ionization rate
    kte = np.arange(1,35,0.01)
    x0 = 7.86
    x1 = 16.37
    x2 = 26.0
    x3 = 38.2
    x4 = 51.6
    invtau0 = 1/(np.exp(x0/kte)/np.sqrt(kte))
    invtau1 = 1/(np.exp(x1/kte)/np.sqrt(kte))
    invtau2 = 1/(np.exp(x2/kte)/np.sqrt(kte))
    invtau3 = 1/(np.exp(x3/kte)/np.sqrt(kte))
    invtau4 = 1/(np.exp(x4/kte)/np.sqrt(kte))
    plt.close()
    plt.plot(kte, invtau0, 'red', label='0 to 1')
    plt.plot(kte, invtau1, 'darkorange', label='1 to 2')
    plt.plot(kte, invtau2, 'gold', label='2 to 3')
    plt.plot(kte, invtau3, 'green', label='3 to 4')
    plt.plot(kte, invtau4, 'blue', label='4 to 5')
    plt.xlabel('Te [eV]')
    plt.ylabel('Relative Ionization Rate [1/s]')
    plt.title('Te Dependence of Ionization Rate for W')
    plt.legend()
    plt.savefig('plots/surface_profiles/IonizRate(Te)')
    
    #ionization rate along surface
    surf_IonizRate0 = ne/(np.exp(x0/te)/np.sqrt(te))
    surf_IonizRate1 = ne/(np.exp(x1/te)/np.sqrt(te))
    surf_IonizRate2 = ne/(np.exp(x2/te)/np.sqrt(te))
    surf_IonizRate3 = ne/(np.exp(x3/te)/np.sqrt(te))
    surf_IonizRate4 = ne/(np.exp(x4/te)/np.sqrt(te))
    plt.close()
    plt.plot(rmrs, surf_IonizRate0, 'red', label='0 to 1')
    plt.plot(rmrs, surf_IonizRate1, 'darkorange', label='1 to 2')
    plt.plot(rmrs, surf_IonizRate2, 'gold', label='2 to 3')
    plt.plot(rmrs, surf_IonizRate3, 'green', label='3 to 4')
    plt.plot(rmrs, surf_IonizRate4, 'blue', label='4 to 5')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Relative Ionization Rate [1/s]')
    plt.title('Ionization Rate along SAS-VW')
    plt.legend()
    plt.savefig('plots/surface_profiles/IonizRate')

    return

def plot_history2D(basic=1, continuousChargeState=0, endChargeState=0):
    history = netCDF4.Dataset("history.nc", "r", format="NETCDF4")

    nP = len(history.dimensions['nP'])
    nT = len(history.dimensions['nT'])
    x = history.variables['x']
    z = history.variables['z']
    charge = history.variables['charge']

    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Particle Tracks')
        
    #define charge state to color mapping
    colors = {0:'black', 1:'red', 2:'orange', 3:'olive', 4:'green', 5:'cyan', \
              6:'purple', 7:'darkmagenta', 8:'pink', 9:'deep pink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:])
    
    if continuousChargeState==1:
        for p in range(15,16):
            print(p,'out of', nP)
            for t in range(0,nT-1):
                plt.plot(x[p][t:t+2],z[p][t:t+2], colors[charge[p][t]])

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:], colors[charge[p][-1]])
        
    plt.xlim(1.43, 1.56)
    plt.ylim(1.07, 1.25)
    plt.savefig('plots/history.pdf')

    return

def plot_surf_nc(pps_per_nP, gitr_rz='../setup/assets/gitr_rz.txt', W_fine_file='../setup/assets/W_fine.txt'):
    surface = netCDF4.Dataset("surface_nP1e5_nT1e6.nc", "r", format="NETCDF4")
    #print(surface.variables['grossErosion'][:])
    
    #calculate area from wall
    #import wall geometry to plot over
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
    
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
    
    area = 0.12*np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))#*np.pi*(r1+r2)

    grossEro = (surface.variables['grossErosion'][:])
    grossDep = (surface.variables['grossDeposition'][:])
    netEro = (grossEro-grossDep)
    
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net eroded comp particles',sum(netEro))
    
    grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    netEro = np.average([netEro[:-1], netEro[1:]],axis=0)/area
    
    print('rmrs length',len(rmrs))
    print('surf length',len(grossEro))
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net eroded flux',sum(netEro))
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netEro_cumsum = np.cumsum(netEro)
    
    plt.close()
    plt.plot(rmrs,grossEro,'r', label='Gross Erosion')
    plt.plot(rmrs,grossDep,'g', label='Redeposition')
    plt.plot(rmrs,netEro,'k', label='Net Erosion')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper left')
    plt.title('GITR Predicted Erosion and Redeposition Profiles')
    plt.savefig('plots/surface.pdf')
    
    plt.close()
    plt.plot(rmrs,grossEro_cumsum,'r', label='Gross Erosion')
    plt.plot(rmrs,grossDep_cumsum,'g', label='Redeposition')
    plt.plot(rmrs,netEro_cumsum,'k', label='Net Erosion')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper left')
    plt.title('GITR Predicted Cumulative Sum of\nErosion and Redeposition Profiles')
    plt.savefig('plots/surface_cumsum.pdf')

def plot_particle_source():
    particleSource = netCDF4.Dataset("../input/particleSource.nc", "r", format="NETCDF4")
    
    x = particleSource.variables['x'][:]
    y = particleSource.variables['y'][:]
    z = particleSource.variables['z'][:]
    
    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.scatter(x,z,marker='_',s=8)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Spatial Particle Source \n nP='+str(len(x)))
    plt.savefig('plots/particleSource.png')


if __name__ == "__main__":
    #init()
    #plot_gitr_gridspace()
    #plot_surf_plasma_params()
    #plot_history2D()
    #plot_surf_nc(1,12,111879178639.80714)
    plot_particle_source()
