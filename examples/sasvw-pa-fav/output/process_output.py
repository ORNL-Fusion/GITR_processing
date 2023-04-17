import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import solps

def init(W_indices = np.arange(10,24)):
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    
    #check which target is the region of interest
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

def plot_history2D(history_file, basic=1, continuousChargeState=0, endChargeState=0):
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")

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
        plt.title('Particle Tracks by End Charge State')
        
    #plt.xlim(1.43, 1.56)
    #plt.ylim(1.07, 1.25)
    plt.savefig('plots/history.pdf')
    plt.close()

    return

def plot_surf_nc(pps_per_nP, \
                 surface_file="surface.nc", \
                 gitr_rz='../setup/assets/gitr_rz.txt', \
                 W_fine_file='../setup/assets/W_fine.txt', \
                 rmrs_fine_file='../setup/assets/rmrs_fine.txt', \
                 plot_cumsum=0):
    
    surface = netCDF4.Dataset(surface_file, "r", format="NETCDF4")
    #print(surface.variables['grossErosion'][:])
    
    #calculate area from wall
    #import wall geometry to plot over
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
    
    #import W surface indices
    with open(W_fine_file, 'r') as file:
        W_fine = file.readlines()
    W_fine = np.array(W_fine,dtype='int')
    
    #import refined rmrs at the W surface
    with open(rmrs_fine_file, 'r') as file:
        rmrs_fine = file.readlines()   
    rmrsFine = np.array(rmrs_fine,dtype='float')
    
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
    
    dist = np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))
    area = np.pi*(r1+r2)*dist

    grossEro = (surface.variables['grossErosion'][:])
    print(grossEro)
    grossDep = (surface.variables['grossDeposition'][:])
    print(grossDep)
    netEro = (grossEro-grossDep)
    
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net eroded comp particles',sum(netEro))
    print('ghost cells',grossEro[-1],grossDep[-1],netEro[-1])
    
    #grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    #grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    grossEro = grossEro[:-1]*pps_per_nP/area
    grossDep = grossDep[:-1]*pps_per_nP/area
    netEro = netEro[:-1]*pps_per_nP/area
    
    print('rmrs length',len(rmrsFine))
    print('surf length',len(grossEro))
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net eroded flux',sum(netEro))
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netEro_cumsum = np.cumsum(netEro)
    
    plt.close()
    plt.plot(rmrsFine,grossEro,'r', label='Gross Erosion')
    plt.plot(rmrsFine,grossDep,'g', label='Redeposition')
    plt.plot(rmrsFine,netEro,'k', label='Net Erosion')
    plt.plot(rmrsFine,np.zeros(len(rmrsFine)),'gray')
    plt.axvline(x=rmrs[4], color='k', linestyle='dotted')
    plt.axvline(x=rmrs[10], color='k', linestyle='dotted')
    plt.axvline(x=rmrs[11], color='k', linestyle='dotted')
    plt.axvline(x=rmrs[12], color='k', linestyle='dotted')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend()#loc='upper left')
    plt.title('GITR Predicted Erosion and \n Redeposition Profiles, nP=1e5, nT=1e6')
    plt.savefig('plots/surface.png')
    
    if plot_cumsum:
        plt.close()
        plt.plot(rmrsFine,grossEro_cumsum,'r', label='Gross Erosion')
        plt.plot(rmrsFine,grossDep_cumsum,'g', label='Redeposition')
        plt.plot(rmrsFine,netEro_cumsum,'k', label='Net Erosion')
        plt.yscale('log')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Flux [#/m2s]')
        plt.legend(loc='upper left')
        plt.title('GITR Predicted Cumulative Sum of \n Erosion and Redeposition Profiles')
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
    
    plt.close()
    plt.hist(z,bins=10)
    plt.xlabel('z [m]')
    plt.title('Spatial Distribution in Z \n nP='+str(len(x)))
    plt.savefig('plots/zhist.png')


if __name__ == "__main__":
    #init()
    #plot_gitr_gridspace()
    plot_history2D("history.nc")
    #plot_surf_nc(37914807680566.16, "surfaces/nP5/surf-5-6.nc")
    #plot_particle_source()
