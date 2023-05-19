import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

def plot_gitr_gridspace():
    profiles, W_indices, R, Z, rmrs = init()
    
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

def plot_particle_source():
    profiles, W_indices, R, Z, rmrs = init()
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

def plot_history2D(history_file, basic=0, continuousChargeState=1, endChargeState=0):
    profiles, W_indices, R, Z, rmrs = init()
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")

    plt.rcParams.update({'lines.linewidth':0.2})
    plt.rcParams.update({'lines.markersize':0})
    plt.rcParams.update({'font.size':16})

    nP = len(history.dimensions['nP'])
    nT = len(history.dimensions['nT'])
    x = history.variables['x']
    z = history.variables['z']
    charge = history.variables['charge']

    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')
    plt.title('W Impurity Trajectories', fontsize=20)
        
    #define charge state to color mapping
    colors = {0:'black', 1:'firebrick', 2:'darkorange', 3:'gold', 4:'limegreen', 5:'dodgerblue', \
              6:'mediumpurple', 7:'darkviolet', 8:'darkmagenta', 9:'deep pink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:])
    
    if continuousChargeState==1:
        for p in np.arange(0,nP,5):
            print(p,'out of', nP)
            for t in np.arange(0,nT,100):
                plt.plot(x[p][t:t+100],z[p][t:t+100], colors[charge[p][t]])

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:], colors[charge[p][-1]])
        plt.title('Particle Tracks by End Charge State')
    
    legend_dict = {'+0':'black', '+1':'firebrick', '+2':'darkorange', '+3':'gold', '+4':'limegreen', '+5':'dodgerblue', \
              '+6':'mediumpurple', '+7':'darkviolet', '+8':'darkmagenta', '+9':'deeppink', '+10':'gray'}
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    plt.legend(handles=patchList, fontsize=12)
    
    plt.xlim(1.3, 1.6)
    plt.ylim(1, 1.3)
    plt.show(block=True)
    plt.savefig('plots/history.pdf')
    plt.close()

    return

def plot_surf_nc(pps_per_nP, \
                 surface_file="surface.nc", \
                 gitr_rz='../setup/assets/gitr_rz.txt', \
                 W_fine_file='../setup/assets/W_fine.txt', \
                 rmrs_fine_file='../setup/assets/rmrs_fine.txt', \
                 plot_cumsum=0):
    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init()
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
    netDep = (grossDep-grossEro)
    
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net deposited comp particles',sum(netDep))
    print('ghost cells',grossEro[-1],grossDep[-1],netDep[-1])
    
    #grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    #grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    grossEro = grossEro[:-1]*pps_per_nP/area
    grossDep = grossDep[:-1]*pps_per_nP/area
    netDep = netDep[:-1]*pps_per_nP/area
    
    print('rmrs length',len(rmrsFine))
    print('surf length',len(grossEro))
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net deposited flux',sum(netDep))
    
    fluxC = 2.22351742230795e22
    grossEro_norm = grossEro/fluxC
    grossDep_norm = grossDep/fluxC
    netDep_norm = netDep/fluxC
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netDep_cumsum = np.cumsum(netDep)
    
    #take gross erosion of a slice in the filterscope range
    r_sp, z_sp = 1.49814916, 1.19640505
    
    r1_start, z1_start = 1.491288, 1.21629
    r1_end, z1_end = 1.49274, 1.22149
    rmrs1_start = -1*np.sqrt((z1_start-z_sp)**2 + (r1_start-r_sp)**2)
    rmrs1_end = -1*np.sqrt((z1_end-z_sp)**2 + (r1_end-r_sp)**2)
    
    r3_start, z3_start = 1.493556, 1.16204
    r3_end, z3_end = 1.490536, 1.1552 
    rmrs3_start = np.sqrt((z3_start-z_sp)**2 + (r3_start-r_sp)**2)
    rmrs3_end = np.sqrt((z3_end-z_sp)**2 + (r3_end-r_sp)**2)
    
    V1_grossEro = 0
    V1_area = 0
    V3_grossEro = 0
    V3_area = 0
    for i,v in enumerate(rmrsFine):
        if v >= rmrs1_end and v <= rmrs1_start:
            V1_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
            V1_area += rmrsFine[i+1] - rmrsFine[i]
        elif v >= rmrs3_start and v <= rmrs3_end:
            V3_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
            V3_area += rmrsFine[i+1] - rmrsFine[i]
    V1_grossEro = V1_grossEro / V1_area
    V3_grossEro = V3_grossEro / V3_area
    
    print('gross erosion in View 1:', V1_grossEro)
    print('gross erosion in View 3:', V3_grossEro)
    
    plt.rcParams.update({'font.size':16})
    plt.rcParams.update({'lines.linewidth':3}) 
    
    plt.close()
    plt.plot(rmrsFine,np.zeros(len(rmrsFine)),'gray')
    plt.axvline(x=rmrs[4], color='k', linestyle='dotted', label='\u0394\u03A8$_B$')
    plt.axvline(x=rmrs[10], color='k', linestyle='dotted')
    plt.axvline(x=rmrs[11], color='k', linestyle='dotted')
    plt.axvline(x=rmrs[12], color='k', linestyle='dotted')
    plt.axvline(x=rmrs1_start, color='lightsalmon', linestyle='dashed', label='Spec View 1')
    plt.axvline(x=rmrs1_end, color='lightsalmon', linestyle='dashed')
    plt.axvline(x=rmrs3_start, color='chocolate', linestyle='dashed', label='Spec View 3')
    plt.axvline(x=rmrs3_end, color='chocolate', linestyle='dashed')
    
    plt.plot(rmrsFine,grossEro_norm,'r', label='Gross Erosion')
    plt.plot(rmrsFine,grossDep_norm,'g', label='Redeposition')
    plt.plot(rmrsFine,netDep_norm,'k', label='Net Deposition')
    
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('\u0393$_{W,outgoing}$ / \u0393$_{C,incoming}$')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    plt.legend(fontsize=10)#loc='upper left')
    plt.title('GITR Predicted Erosion and \n Redeposition Profiles, nP=1e6, nT=1e6')
    plt.savefig('plots/surface.png')
    
    if plot_cumsum:
        plt.close()
        plt.plot(rmrsFine,grossEro_cumsum,'r', label='Gross Erosion')
        plt.plot(rmrsFine,grossDep_cumsum,'g', label='Redeposition')
        plt.plot(rmrsFine,netDep_cumsum,'k', label='Net Erosion')
        plt.yscale('log')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Flux [m$^{-2}$s$^{-1}$]')
        plt.legend(loc='upper left')
        plt.title('GITR Predicted Cumulative Sum of \n Erosion and Redeposition Profiles',fontsize=20)
        plt.savefig('plots/surface_cumsum.pdf')
    

def spectroscopy(pps_per_nP, View=3, \
                 specFile='spec.nc'):
    
    spec = netCDF4.Dataset(specFile, "r", format="NETCDF4")
    profiles, W_indices, R, Z, rmrs = init()
    
    gridr = spec.variables['gridR'][:]
    gridz = spec.variables['gridZ'][:]
    dr = gridr[1]-gridr[0]
    dz = gridz[1]-gridz[0]
    gridr = np.append(gridr, gridr[-1]+dr)
    gridz = np.append(gridz, gridz[-1]+dz)
    rr, zz = np.meshgrid(gridr,gridz)
    density_unitless = spec.variables['n'][:]
    density_neutrals = density_unitless[0]
    #density_neutrals = np.sum(density_unitless, axis=0)
    
    plt.pcolor(gridr,gridz,density_neutrals)
    plt.axis('Scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar(label='# Computational Particles')
    plt.title('Raw W0 Spec Output form GITR')
    plt.savefig('plots/spec_compParticles.png')
    
    filterscope_diameter = 0.00702068
    tokamak_circumference = 2 * np.pi * 1.4925
    toroidal_fraction = filterscope_diameter / tokamak_circumference
    
    r_mid = dr/2 + rr[:-1,:]
    dV = 2 * np.pi * r_mid[:,:-1] * dr * dz * toroidal_fraction
    # dA = np.pi * (0.00702068/2)**2
    density_volumetric = pps_per_nP * (density_neutrals/dV) * 1e-8 * 3e3 # W0 m-2 s-1
    #neutral_flux = density_volumetric * 3e3 # W0 m-2 s-1
    #print('TESTTT',neutral_flux)
    
    plt.close()
    plt.plot(R,Z)
    plt.pcolor(gridr,gridz,density_volumetric)
    plt.axis('Scaled')
    plt.xlim(min(R),max(R))
    plt.ylim(min(Z),max(Z))
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
    plt.title('Toroidally Integrated W0 Density')
    plt.savefig('plots/spec_density.png')

    view_fraction = toroidal_fraction * np.pi/4
    density_sliced = view_fraction * density_volumetric
    
    if View==3:
        line1 = -0.97*gridr + 2.600823
        line2 = -0.906*gridr + 2.515464
        
        if True:
            plt.close()
            plt.plot(R,Z)
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,density_sliced)
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('r [m]')
            plt.ylabel('z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
            plt.title('Toroidal Slice of W0 Density')
            plt.savefig('plots/spec_density_sliced.png')
        
        rstart, rend = 39, 72
        zstart, zend = 40, 80
        gridr = gridr[rstart:rend]
        gridz = gridz[zstart:zend]
        fscope = density_sliced[zstart:zend, rstart:rend]
        
        fscope[np.where(fscope==0)] = 'NaN'
        #fscope[-5:,23:] = 'NaN'
    
        line1 = -0.97*gridr + 2.600823
        line2 = -0.906*gridr + 2.515464
        
        if True:
            plt.close()
            plt.plot(R,Z)
            line1 = -0.97*gridr + 2.600823
            line2 = -0.906*gridr + 2.515464    
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,fscope,shading='auto')
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('r [m]')
            plt.ylabel('z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
    
    fscope[np.isnan(fscope)] = 0
    print('W0 Neutral Flux:', np.sum(fscope)) # in units of m-3 s-1
    plt.title('Toroidal Slice of W0 Density \n W0: %10e' %np.sum(fscope))
    plt.savefig('plots/spec_filterscope.png')
    
if __name__ == "__main__":
    #init()
    #plot_gitr_gridspace()
    #plot_particle_source()
    #plot_history2D("../../../../GITR/scratch/output/history.nc")
    plot_surf_nc(3791480768056.615, "surfaceP6T6.nc")
    #plot_surf_nc(37914807680566.16, "/Users/Alyssa/Dev/SAS-VW-Data/netcdf_data/nP5/surf-5-6.nc")
    #spectroscopy(3791480768056.615,specFile='specP6T6.nc')
