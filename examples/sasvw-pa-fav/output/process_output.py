import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))
sys.path.insert(0, os.path.abspath('../setup/'))

import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as path
import netCDF4
import solps
import makeParticleSource

def init(W_indices = np.arange(11,22), plot_rz=False):
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    
    #check which target is the region of interest
    r_inner_target = profiles.variables['r_inner_target'][:]
    z_inner_target = profiles.variables['z_inner_target'][:]
    rmrs = profiles.variables['rmrs_inner_target_midpoints'][W_indices]
    
    #plot r,z to check that target indices are restricted to W surfaces
    if plot_rz:
        plt.close()
        plt.plot(r_inner_target, z_inner_target)
    
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':2})
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
    plt.hist(z,bins=10)
    plt.xlabel('z [m]')
    plt.title('Spatial Distribution in Z \n nP='+str(len(x)))
    plt.savefig('plots/zhist.png')
    
    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.scatter(x,z,marker='_')
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Spatial Particle Source \n nP='+str(len(x)))
    plt.savefig('plots/particleSource.png')

def plot_history2D(history_file='history.nc', \
                   basic=0, continuousChargeState=1, endChargeState=0, \
                   plot_particle_source=0, markersize=0):
    
    if plot_particle_source:
        particleSource = netCDF4.Dataset("../input/particleSource.nc", "r", format="NETCDF4")
        x0 = particleSource.variables['x'][:]
        z0 = particleSource.variables['z'][:]
    
    profiles, W_indices, R, Z, rmrs = init()
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")

    plt.rcParams.update({'lines.linewidth':0.3})
    plt.rcParams.update({'lines.markersize':markersize})
    plt.rcParams.update({'font.size':16})

    nP = len(history.dimensions['nP'])
    print('nP:',nP)
    nT = len(history.dimensions['nT'])
    x = history.variables['x'][:]
    y = history.variables['y'][:]
    z = history.variables['z'][:]
    r = np.sqrt(x**2 + y**2)
    charge = history.variables['charge']

    #plt.close()
    plt.close()
    if plot_particle_source: plt.scatter(x0,z0,marker='o',s=10)
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('W Impurity Trajectories', fontsize=20)
        
    #define charge state to color mapping
    colors = {0:'black', 1:'firebrick', 2:'darkorange', 3:'gold', 4:'limegreen', 5:'dodgerblue', \
              6:'mediumpurple', 7:'darkviolet', 8:'darkmagenta', 9:'deep pink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            print(p,'out of', nP)
            plt.plot(r[p][:],z[p][:])
            #plt.scatter(x[p][:],z[p][:],marker='_',s=50,c='k')
    
    if continuousChargeState==1:
        for p in np.arange(0,nP,10):
            print(p,'out of', nP)
            t=0
            while t<nT-1:
                if r[p][t] == r[p][t+1]: t=nT
                else:
                    plt.plot(r[p][t:t+2],z[p][t:t+2], colors[charge[p][t]])
                    t+=1

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(r[p][:],z[p][:], colors[charge[p][-1]])
        plt.title('Particle Tracks by End Charge State')
    
    plt.scatter(1.49829829, 1.19672716, label='Strikepoint', marker='X', color='k', s=100, zorder=5)
    
    legend_dict = {'+0':'black', '+1':'firebrick', '+2':'darkorange', '+3':'gold', '+4':'limegreen', '+5':'dodgerblue', \
              '+6':'mediumpurple', '+7':'darkviolet', '+8':'darkmagenta', '+9':'deeppink', '+10':'gray'}
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    if basic==0: plt.legend(handles=patchList, fontsize=12)
    
    plt.xlim(1.45, 1.525)
    plt.ylim(1.1, 1.23)
    plt.show(block=True)
    plt.savefig('plots/history.pdf')
    plt.close()

    return

def plot_surf_nc(pps_per_nP, nP10, nT10, \
                 tile_shift_indices = [], Bangle_shift_indices = [], \
                 surface_file="surface.nc", \
                 gitr_rz='../setup/assets/gitr_rz.txt', \
                 W_fine_file='../setup/assets/W_fine.txt', \
                 rmrs_fine_file='../setup/assets/rmrs_fine.txt', \
                 norm=None, plot_cumsum=0):
    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init()
    rmrsCoords = profiles.variables['rmrs_inner_target'][W_indices]
    surface = netCDF4.Dataset(surface_file, "r", format="NETCDF4")
    partSource_flux, fluxD, fluxC = makeParticleSource.distributed_source(nP=int(10**nP10), \
                surfW=np.arange(11,22), \
                tile_shift_indices = [1,9], \
                Bangle_shift_indices = [3,8,9])
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
    
    print('\n')
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net deposited comp particles',sum(netDep))
    print('ghost cells',grossEro[-1],grossDep[-1],netDep[-1])
    
    #grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    #grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    grossEro = grossEro[:-1]*pps_per_nP/area
    grossDep = grossDep[:-1]*pps_per_nP/area
    netDep = netDep[:-1]*pps_per_nP/area
    
    print('\n')
    print('rmrs length',len(rmrsFine))
    print('surf length',len(grossEro))
    print('\n')
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net deposited flux',sum(netDep))
    print('self-sputtering fraction',(sum(grossEro)-sum(partSource_flux))/sum(grossEro))
    
    if norm=='C':
        grossEro_norm = grossEro/fluxC
        grossDep_norm = grossDep/fluxC
        netDep_norm = netDep/fluxC
    elif norm=='D':
        grossEro_norm = grossEro/fluxD
        grossDep_norm = grossDep/fluxD
        netDep_norm = netDep/fluxD
    else: 
        grossEro_norm = grossEro
        grossDep_norm = grossDep
        netDep_norm = netDep
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netDep_cumsum = np.cumsum(netDep)
    
    #take gross erosion of a slice in the filterscope range
    r_sp, z_sp = 1.49829829, 1.19672716
    
    r1_start, z1_start = 1.491288, 1.21629
    r1_end, z1_end = 1.49274, 1.22149
    rmrs1_start = -1*np.sqrt((z1_start-z_sp)**2 + (r1_start-r_sp)**2)
    rmrs1_end = -1*np.sqrt((z1_end-z_sp)**2 + (r1_end-r_sp)**2)
    
    r2_start, z2_start = 1.492041, 1.19418 
    r2_end, z2_end = 1.492716, 1.18709 
    rmrs2_start = 1*np.sqrt((z2_start-z_sp)**2 + (r2_start-r_sp)**2)
    rmrs2_end = 1*np.sqrt((z2_end-z_sp)**2 + (r2_end-r_sp)**2)
    
    r3_start, z3_start = 1.493556, 1.16204
    r3_end, z3_end = 1.490536, 1.1552 
    rmrs3_start = np.sqrt((z3_start-z_sp)**2 + (r3_start-r_sp)**2)
    rmrs3_end = np.sqrt((z3_end-z_sp)**2 + (r3_end-r_sp)**2)
    
    #initialize intermediate variables for calculating erosion within a given view
    V1_indices = np.empty(0, dtype=int)
    V2_indices = np.empty(0, dtype=int)
    V3_indices = np.empty(0, dtype=int)
    V1_grossEro = 0
    V2_grossEro = 0
    V3_grossEro = 0
    V1_area = 0
    V2_area = 0
    V3_area = 0
    
    for i,v in enumerate(rmrsFine):
        if v >= rmrs1_end and v <= rmrs1_start:
            V1_indices = np.append(V1_indices, i)
        elif v >= rmrs2_start and v <= rmrs2_end:
            V2_indices = np.append(V2_indices, i)
        elif v >= rmrs3_start and v <= rmrs3_end:
            V3_indices = np.append(V3_indices, i)
    
    #add all cells beginning inside the view
    for i in V1_indices:
        V1_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
        V1_area += rmrsFine[i+1] - rmrsFine[i]
    for i in V2_indices:
        V2_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
        V2_area += rmrsFine[i+1] - rmrsFine[i]
    for i in V3_indices:
        V3_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
        V3_area += rmrsFine[i+1] - rmrsFine[i]
    
    #add fraction of previous cell ending inside the view
    V1_grossEro += grossEro[V1_indices[0]-1] * (rmrs1_end-rmrsFine[V1_indices[0]])
    V1_area += rmrs1_end-rmrsFine[V1_indices[0]]
    V2_grossEro += grossEro[V2_indices[0]-1] * (rmrs2_start-rmrsFine[V2_indices[0]])
    V2_area += rmrs2_start-rmrsFine[V2_indices[0]]
    V3_grossEro += grossEro[V3_indices[0]-1] * (rmrs3_start-rmrsFine[V3_indices[0]])
    V3_area += rmrs3_start-rmrsFine[V3_indices[0]]

    #subtract fraction of last cell starting inside the view
    V1_grossEro -= grossEro[V1_indices[-1]] * (rmrs1_start-rmrsFine[V1_indices[-1]+1])
    V1_area -= rmrs1_start-rmrsFine[V1_indices[-1]+1] 
    V2_grossEro -= grossEro[V2_indices[-1]] * (rmrs2_end-rmrsFine[V2_indices[-1]+1])
    V2_area -= rmrs2_end-rmrsFine[V2_indices[-1]+1]
    V3_grossEro -= grossEro[V3_indices[-1]] * (rmrs3_end-rmrsFine[V3_indices[-1]+1])
    V3_area -= rmrs3_end-rmrsFine[V3_indices[-1]+1]
    
    V1_grossEro = V1_grossEro / V1_area
    V2_grossEro = V2_grossEro / V2_area
    V3_grossEro = V3_grossEro / V3_area
    
    print('\n')
    print('gross erosion in View 1:', V1_grossEro)
    print('gross erosion in View 2:', V2_grossEro)
    print('gross erosion in View 3:', V3_grossEro)
    
    plt.rcParams.update({'font.size':30})
    plt.rcParams.update({'lines.linewidth':5}) 
    
    #plot self-sputtering
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsFine, 100*(grossEro-partSource_flux)/grossEro)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Percentage')
    plt.title('Percentage of Gross Erosion from Self-Sputtering')
    plt.show(block=True)
    
    #plot main surface plot with 3 views
    plt.close()
    plt.axvspan(rmrs1_start, rmrs1_end, color='#f7bc00', alpha=0.5)
    plt.axvspan(rmrs2_start, rmrs2_end, color='lightsalmon', alpha=0.5)
    plt.axvspan(rmrs3_start, rmrs3_end, color='#f99301', alpha=0.5)
    
    plt.plot(rmrsFine,np.zeros(len(rmrsFine)),'gray')
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsFine,grossEro_norm,'r', label='Gross Erosion')
    plt.plot(rmrsFine,grossDep_norm,'g', label='Redeposition')
    plt.plot(rmrsFine,netDep_norm,'k', label='Net Deposition')
    
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('\u0393$_{W,outgoing}$ / \u0393$_{C,incoming}$')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    plt.legend(fontsize=20)#loc='upper left')
    plt.title('GITR Predicted Erosion and \n Redeposition Profiles, nP=1e'+str(nP10)+', nT=1e'+str(nT10), fontsize=30)
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
    
    return
        
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
    
def analyze_leakage(historyFile):
    bFile = '../input/bField.nc'
    bField = netCDF4.Dataset(bFile)
    r_bField = bField.variables['r'][:]
    z_bField = bField.variables['z'][:]
    psi = bField.variables['psi'][:]
    xpoint_z = 0.9874948
    
    gitr_rz='../setup/assets/gitr_rz.txt'
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
        
    r_wall = np.zeros(len(wall))
    z_wall = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        r_wall[i] = float(point[0])
        z_wall[i] = float(point[1])
    
    history = netCDF4.Dataset(historyFile)
    nP = len(history.dimensions['nP'])
    x = history.variables['x'][:]
    y = history.variables['y'][:]
    z = history.variables['z'][:]
    r = np.sqrt(x**2 + y**2)
    
    cs = plt.contourf(r_bField,z_bField,psi,1)
    p = cs.collections[0].get_paths()[0]
    v = p.vertices
    pathx = v[:,0]
    pathy = v[:,1]
    
    plt.close()
    fig,ax = plt.subplots()
    ax.fill(pathx,pathy,'lightpink')
    plt.plot(r_wall, z_wall, 'darkmagenta', linewidth=2)
    for p in range(0,nP,100):
        plt.plot(r[p][:],z[p][:],'k',linewidth=0.3)
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Leakage of W into the core')
    plt.axis('Scaled')
    
    leakage = 0
    polygon = path.Path(v,closed=True)
    for i in range(0,nP):
        if polygon.contains_point((r[i][-1],z[i][-1])) and z[i][-1]<=xpoint_z: leakage+=1
    
    print('leakage fraction:', leakage/nP)
    
    
    
    return    

def analyze_forces(varString, component, rzlim=True, colorbarLimits=[], dt=1e-8):
    #import wall geometry to plot over
    gitr_rz='../setup/assets/gitr_rz.txt'
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
        
    r_wall = np.zeros(len(wall))
    z_wall = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        r_wall[i] = float(point[0])
        z_wall[i] = float(point[1])
    
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    gridr = profiles.variables['gridr'][:]
    gridz = profiles.variables['gridz'][:]
    
    Br = profiles.variables['Br'][:]
    Bt = profiles.variables['Bt'][:]
    Bz = profiles.variables['Bz'][:]
    Er = profiles.variables['Er'][:]
    Et = profiles.variables['Et'][:]
    Ez = profiles.variables['Ez'][:]
    charge = profiles.variables['charge'][:]
    density = profiles.variables['ni'][:]
    ti = profiles.variables['ti'][:]
    te = profiles.variables['te'][:]
    vr_background = profiles.variables['vr'][:]
    vt_background = profiles.variables['vt'][:]
    vz_background = profiles.variables['vz'][:]
    
    amu = 183.84
    amu_D = 2
    amu_C = 12.011
    Q = 1.60217662e-19
    
    charge_W = np.ones(np.shape(Br))
    q = Q*charge_W
    vr = -1000*np.ones(np.shape(Br))
    vt = -500*np.ones(np.shape(Br))
    vz = 700*np.ones(np.shape(Br))
    
    if varString=='q(v x B)':
        vartype = 'F'
        titleString = ' = (q(v x B))'
        
        Fr = q*(vt*Bz - vz*Bt)
        Ft = q*(vz*Br - vr*Bz)
        Fz = q*(vr*Bt - vt*Br)

    
    elif varString=='Lorentz' or varString=='Lorentz dv':
        vartype = 'F'
        titleString = ' = (q(E + v x B))'
        
        Fr = q*(Er + vt*Bz - vz*Bt)
        Ft = q*(Et + vz*Br - vr*Bz)
        Fz = q*(Ez + vr*Bt - vt*Br)
    
        if varString=='Lorentz dv':
            vartype = 'v'
            Wmass_amu = 183.84
            Wmass_kg = Wmass_amu * 1.66054e-27
            
            Fr = dt*Fr/Wmass_kg
            Ft = dt*Ft/Wmass_kg
            Fz = dt*Fz/Wmass_kg
    
    elif varString=='ExB drift':
        vartype = 'v'
        titleString = ' = ($v_{ExB}$)'
        
        dr = np.average(gridr[1:]-gridr[:-1])
        dz = np.average(gridz[1:]-gridz[:-1])
        Wmass_amu = 183.84
        Wmass_kg = Wmass_amu * 1.66054e-27
        B_mag2 = Br**2 + Bt**2 + Bz**2
        B_mag = np.sqrt(B_mag2)
        vdotB = vr*Br + vt*Bt + vz*Bz
        v_parallel_r = Br * vdotB / B_mag2
        v_parallel_t = Bt * vdotB / B_mag2
        v_parallel_z = Bz * vdotB / B_mag2
        v_parallel = np.sqrt(v_parallel_r**2 + v_parallel_t**2 + v_parallel_z**2)
        v_mag = np.sqrt(vr**2 + vt**2 + vz**2)
        v_perp = v_mag - v_parallel
        r_Larmor = Wmass_kg * v_perp / (q * B_mag)
        print('r_Larmor:', np.average(r_Larmor))
        
        ExB_r = Et*Bz - Ez*Bt
        ExB_t = Ez*Br - Er*Bz
        ExB_z = Er*Bt - Et*Br
        
        del_ExB_z = np.zeros(np.shape(Br))
        del_ExB_z[0,:] = (ExB_z[1,:] - ExB_z[0,:])/dz
        del_ExB_z[-1,:] = (ExB_z[-1,:] - ExB_z[-2,:])/dz
        for i in range(1,len(gridz)-1):
            del_ExB_z[i,:] = (ExB_z[i+1,:] - ExB_z[i-1,:])/(2*dz)
        
        del2_ExB_z = np.zeros(np.shape(Br))
        del2_ExB_z[0,:] = (del_ExB_z[1,:] - del_ExB_z[0,:])/dz
        del2_ExB_z[-1,:] = (del_ExB_z[-1,:] - del_ExB_z[-2,:])/dz
        for i in range(1,len(gridz)-1):
            del2_ExB_z[i,:] = (del_ExB_z[i+1,:] - del_ExB_z[i-1,:])/(2*dz)
        
        del_ExB_r = np.zeros(np.shape(Br))
        del_ExB_r[:,0] = gridr[0] * (ExB_r[:,1] - ExB_r[:,0])/dr
        del_ExB_r[:,-1] = gridr[-1] * (ExB_r[:,-1] - ExB_r[:,-2])/dr
        for i in range(1,len(gridr)-1):
            del_ExB_r[:,i] = gridr[i] * (ExB_r[:,i+1] - ExB_r[:,i-1])/(2*dr)
        
        del2_ExB_r = np.zeros(np.shape(Br))
        del2_ExB_r[:,0] = (del_ExB_r[:,1] - del_ExB_r[:,0])/(gridr[0]*dr)
        del2_ExB_r[:,-1] = (del_ExB_r[:,-1] - del_ExB_r[:,-2])/(gridr[-1]*dr)
        for i in range(1,len(gridr)-1):
            del2_ExB_r[:,i] = (del_ExB_r[:,i+1] - del_ExB_r[:,i-1])/(2*gridr[i]*dr)
        
        Fr = ExB_r/B_mag2 + (r_Larmor**2 * del2_ExB_r)/(4 * B_mag2)
        Ft = ExB_t/B_mag2
        Fz = ExB_z/B_mag2 + (r_Larmor**2 * del2_ExB_z)/(4 * B_mag2)
    
    elif varString=='drag' or varString=='drag dv':
        vartype = 'F'
        titleString = ' = ($m \\nu_S U_{\parallel}$)'
        
        Z_D = 1
        Z_C = 6
        EPS0 = 8.854187e-12 #epsilon_0 = vacuum permissivity
        MI = 1.6737236e-27
        
        #inverese acceleration of D and C background ions
        a_D = amu_D*MI/(2*ti*Q)
        a_C = amu_C*MI/(2*ti*Q)
        
        relativeVelocity_r = vr - vr_background
        relativeVelocity_t = vt - vt_background
        relativeVelocity_z = vz - vz_background
        velocityNorm = np.sqrt(relativeVelocity_r**2 + relativeVelocity_t**2 + relativeVelocity_z**2)
        
        #plasma parameter
        lam_d_D = np.sqrt(EPS0*te/(density*(Z_D**2)*Q))
        lam_d_C = np.sqrt(EPS0*te/(density*(Z_C**2)*Q))
        lam_D = 12.0*np.pi*density*lam_d_D**3/charge
        lam_C = 12.0*np.pi*density*lam_d_C**3/charge
        gam_electron_background_D = 0.238762895*charge**2*np.log(lam_D)/(amu*amu)
        gam_electron_background_C = 0.238762895*charge**2*np.log(lam_C)/(amu*amu)
        nu_0_D = gam_electron_background_D*density/velocityNorm**3
        nu_0_C = gam_electron_background_C*density/velocityNorm**3
        
        #wtf is any of this
        xx_D = velocityNorm**2 * a_D
        xx_C = velocityNorm**2 * a_C
        psi_prime_D = 2.0*np.sqrt(xx_D/np.pi)*np.exp(-xx_D)
        psi_prime_C = 2.0*np.sqrt(xx_C/np.pi)*np.exp(-xx_C)
        psi_psiprime_D = special.erf(np.sqrt(xx_D))
        psi_psiprime_C = special.erf(np.sqrt(xx_C))
        psi_D = psi_psiprime_D - psi_prime_D
        psi_C = psi_psiprime_C - psi_prime_C
        
        #Find zero velocity moment friction frequency:
        #nu_0 for W in D (nu_friction_D) and for W in C (nu_friction_C)
        nu_friction_D = (1+amu/amu_D)*psi_D*nu_0_D
        nu_friction_C = (1+amu/amu_C)*psi_C*nu_0_C
        nu_s = 0.9*nu_friction_D+0.1*nu_friction_C
        
        B_mag2 = Br**2 + Bt**2 + Bz**2
        vdotB = relativeVelocity_r*Br + relativeVelocity_t*Bt + relativeVelocity_z*Bz
        v_parallel_r = Br * vdotB / B_mag2
        v_parallel_t = Bt * vdotB / B_mag2
        v_parallel_z = Bz * vdotB / B_mag2
        
        #Calulate drag force component of the Fokker-Plank collisional forces
        Fr = -1 * amu * MI * nu_s * v_parallel_r
        Ft = -1 * amu * MI * nu_s * v_parallel_t
        Fz = -1 * amu * MI * nu_s * v_parallel_z
        
        if varString=='drag dv':
            vartype = 'v'
            Wmass_amu = 183.84
            Wmass_kg = Wmass_amu * 1.66054e-27
            
            Fr = dt*Fr/Wmass_kg
            Ft = dt*Ft/Wmass_kg
            Fz = dt*Fz/Wmass_kg
        
    elif varString=='friction':
        vartype = 'F'
        titleString = ' = ($F_{FP}$)'
        
        Z_D = 1
        Z_C = 6
        EPS0 = 8.854187e-12 #epsilon_0 = vacuum permissivity
        MI = 1.6737236e-27
        
        #inverese acceleration of D and C background ions
        a_D = amu_D*MI/(2*ti*Q)
        a_C = amu_C*MI/(2*ti*Q)
        
        relativeVelocity_r = vr - vr_background
        relativeVelocity_t = vt - vt_background
        relativeVelocity_z = vz - vz_background
        velocityNorm = np.sqrt(relativeVelocity_r**2 + relativeVelocity_t**2 + relativeVelocity_z**2)
        
        #wtf is this
        lam_d_D = np.sqrt(EPS0*te/(density*(Z_D**2)*Q))
        lam_d_C = np.sqrt(EPS0*te/(density*(Z_C**2)*Q))
        lam_D = 12.0*np.pi*density*lam_d_D**3/charge
        lam_C = 12.0*np.pi*density*lam_d_C**3/charge
        gam_electron_background_D = 0.238762895*charge**2*np.log(lam_D)/(amu*amu)
        gam_electron_background_C = 0.238762895*charge**2*np.log(lam_C)/(amu*amu)
        nu_0_D = gam_electron_background_D*density/velocityNorm**3
        nu_0_C = gam_electron_background_C*density/velocityNorm**3
        
        #wtf is any of this
        xx_D = velocityNorm**2 * a_D
        xx_C = velocityNorm**2 * a_C
        psi_prime_D = 2.0*np.sqrt(xx_D/np.pi)*np.exp(-xx_D)
        psi_prime_C = 2.0*np.sqrt(xx_C/np.pi)*np.exp(-xx_C)
        psi_psiprime_D = special.erf(np.sqrt(xx_D))
        psi_psiprime_C = special.erf(np.sqrt(xx_C))
        psi_D = psi_psiprime_D - psi_prime_D
        psi_C = psi_psiprime_C - psi_prime_C
        
        #Find zero velocity moment friction frequency:
        #nu_0 for W in D (nu_friction_D) and for W in C (nu_friction_C)
        nu_friction_D = (1+amu/amu_D)*psi_D*nu_0_D
        nu_friction_C = (1+amu/amu_C)*psi_C*nu_0_C
        nu_s = 0.9*nu_friction_D+0.1*nu_friction_C
        
        B_mag2 = Br**2 + Bt**2 + Bz**2
        vdotB = relativeVelocity_r*Br + relativeVelocity_t*Bt + relativeVelocity_z*Bz
        v_parallel_r = Br * vdotB / B_mag2
        v_parallel_t = Bt * vdotB / B_mag2
        v_parallel_z = Bz * vdotB / B_mag2
        
        #Calulate drag force component of the Fokker-Plank collisional forces
        drag_r = -1 * amu * MI * nu_s * v_parallel_r
        drag_t = -1 * amu * MI * nu_s * v_parallel_t
        drag_z = -1 * amu * MI * nu_s * v_parallel_z
        
        nu_parallel_D = psi_D/xx_D*nu_0_D
        nu_parallel_C = psi_C/xx_C*nu_0_C
        nu_p = 0.9*nu_parallel_D+0.1*nu_parallel_C
        nu_deflection_D = 2*(psi_psiprime_D - psi_D/(2*xx_D))*nu_0_D
        nu_deflection_C = 2*(psi_psiprime_C - psi_C/(2*xx_C))*nu_0_C
        nu_d = 0.9*nu_deflection_D+0.1*nu_deflection_C
        nu_energy_D = 2*(amu/amu_D*psi_D - psi_prime_D)*nu_0_D
        nu_energy_C = 2*(amu/amu_C*psi_C - psi_prime_C)*nu_0_C
        nu_e = 0.9*nu_energy_D+0.1*nu_energy_C
        
        print('s:', np.average(nu_friction_D), np.average(nu_friction_C))
        print('parallel:', np.average(nu_parallel_D), np.average(nu_parallel_C))
        print('deflection:', np.average(nu_deflection_D), np.average(nu_deflection_C))
        print('energy:', np.average(nu_energy_D), np.average(nu_energy_C))
        
        v_perp_r = np.abs(relativeVelocity_r - v_parallel_r)
        v_perp_t = np.abs(relativeVelocity_t - v_parallel_t)
        v_perp_z = np.abs(relativeVelocity_z - v_parallel_z)
        
        m = amu * MI
        Fr = drag_r - m*v_parallel_r*np.sqrt(nu_p/dt) - m*v_perp_r*np.sqrt(nu_d/(2*dt))
        Ft = drag_t - m*v_parallel_t*np.sqrt(nu_p/dt) - m*v_perp_t*np.sqrt(nu_d/(2*dt))
        Fz = drag_z - m*v_parallel_z*np.sqrt(nu_p/dt) - m*v_perp_z*np.sqrt(nu_d/(2*dt))
        
        
    elif varString=='gradT' or varString=='gradT dv':
        vartype = 'F'
        titleString = ' = ($ \\alpha \\nabla_{\parallel} T_e + \\beta \\nabla_{\parallel} T_i $)'
        
        gradTir = profiles.variables['gradTir'][:]*Q
        gradTit = profiles.variables['gradTit'][:]*Q
        gradTiz = profiles.variables['gradTiz'][:]*Q
        gradTer = profiles.variables['gradTer'][:]*Q
        gradTet = profiles.variables['gradTet'][:]*Q
        gradTez = profiles.variables['gradTez'][:]*Q
        
        mu_D = amu / (amu_D + amu)
        mu_C = amu / (amu_C + amu)
        alpha = 0.71*charge_W**2
        beta_D = 3 * (mu_D + 5*np.sqrt(2.0) * charge_W**2 * (1.1*mu_D**(5 / 2) - 0.35*mu_D**(3 / 2)) - 1) / (2.6 - 2*mu_D + 5.4*mu_D**2)
        beta_C = 3 * (mu_C + 5*np.sqrt(2.0) * charge_W**2 * (1.1*mu_C**(5 / 2) - 0.35*mu_C**(3 / 2)) - 1) / (2.6 - 2*mu_C + 5.4*mu_C**2)
        #assume a 3.7% C plasma everywhere
        beta = 0.9565*beta_D + 0.0435*beta_C 
        
        Fr = alpha*gradTer + beta*gradTir
        Ft = alpha*gradTet + beta*gradTit
        Fz = alpha*gradTez + beta*gradTiz
        
        if varString=='gradT dv':
            vartype = 'v'
            Wmass_amu = 183.84
            Wmass_kg = Wmass_amu * 1.66054e-27
            
            Fr = dt*Fr/Wmass_kg
            Ft = dt*Ft/Wmass_kg
            Fz = dt*Fz/Wmass_kg
    
        
    Fr[np.where(Fr==0)] = 'nan'
    Ft[np.where(Ft==0)] = 'nan'
    Fz[np.where(Fz==0)] = 'nan'
    
    gridrz = [gridr, gridz, r_wall, z_wall]
    if component == 'r': plot_forces(Fr, vartype+'r'+titleString+'$_r$', gridrz, vartype, rzlim, colorbarLimits)
    if component == 't': plot_forces(Ft, vartype+'t'+titleString+'$_t$', gridrz, vartype, rzlim, colorbarLimits)
    if component == 'z': plot_forces(Fz, vartype+'z'+titleString+'$_z$', gridrz, vartype, rzlim, colorbarLimits)
    
    return 

def plot_forces(var, titleString, gridrz, vartype='F', rzlim=True, colorbarLimits=[]):
    [gridr, gridz, r_wall, z_wall] = gridrz
    plt.rcParams.update({'pcolor.shading':'auto'})
    plt.rcParams.update({'image.cmap':'rainbow'})
    rlim = [1.375, 1.575]
    zlim = [1.05, 1.25]
    
    if rzlim:
        r_indices = np.empty(0, dtype='int')
        z_indices = np.empty(0, dtype='int')
        
        for i,v in enumerate(gridr):
            if v>=rlim[0] and v<=rlim[1]:
                r_indices = np.append(r_indices, i)
        
        for i,v in enumerate(gridz):
            if v>=zlim[0] and v<=zlim[1]:
                z_indices = np.append(z_indices, i)

        gridr = gridr[r_indices]
        gridz = gridz[z_indices]
        var = var[z_indices[0]:z_indices[-1]+1,r_indices[0]:r_indices[-1]+1]
    
    plt.close()
    plt.plot(r_wall,z_wall,'-k')
    plot = plt.pcolor(gridr,gridz,var)
    
    if colorbarLimits != []: plot.set_clim(vmin=colorbarLimits[0],vmax=colorbarLimits[1])
    if vartype=='F': plt.colorbar(label='\n Force [N]')
    if vartype=='v': plt.colorbar(label='\n Velocity [m/s]')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title(titleString)
    plt.axis('Scaled')
    plt.xlim(rlim)
    plt.ylim(zlim)
    
    return

def theoretical_sheath(profiles, W_surf):
    Bmag = profiles.variables['Bmag_inner_target'][W_surf]
    te = profiles.variables['te_inner_target'][W_surf]
    ti = profiles.variables['ti_inner_target'][W_surf]
    ni = profiles.variables['ni_inner_target'][:,W_surf]
    niC = ni[2] + ni[3] + ni[4] + ni[5] + ni[6] + ni[7] + ni[8]
    ni = ni[0] + ni[1] + niC
    
    #Debye Length
    epsilon_0 = 55.26349406
    lambda_D = np.sqrt(epsilon_0 * te / (1e-18 * ni)) #in microns
    lambda_D = lambda_D*1e-6 #in meters
    
    #Larmor Radius
    unit_charge = 1.602176634e-19
    amu2kg = 1.6605402e-27
    mi = 2 * amu2kg
    cs = np.sqrt(unit_charge*(te+ti)/mi) #sound speed
    r_Larmor = (mi/(unit_charge*Bmag)) * cs * 1000 #in milimeters
    r_Larmor = r_Larmor*1e-3 #in meters
    
    #Sheath Width Factors
    Debye_factor = 3
    Chodura_factor = 3
    Debye_width = Debye_factor * lambda_D
    Chodura_width = Chodura_factor * r_Larmor
    
    return Debye_width, Chodura_width


def ionization_analysis(plotting, historyFile, positionsFile, tile_shift_indices, Bangle_shift_indices, W_surf):    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init(W_surf)
    history = netCDF4.Dataset(historyFile)
    positions = netCDF4.Dataset(positionsFile)
    
    z_inner_target = z_inner_target[W_surf]
    rmrsMid = rmrs[:-1]
    rmrsCoords = profiles.variables['rmrs_inner_target'][W_surf]
    nP = len(history.dimensions['nP'])
    charge = history.variables['charge'][:]
    angle = positions.variables['angle'][:]
    
    Debye, Chodura = theoretical_sheath(profiles,W_surf)
        
    x = history.variables['x'][:]
    y = history.variables['y'][:]
    z = history.variables['z'][:]
    
    tally_never_ionizes = np.zeros(len(z_inner_target)-1)
    avg_distance_to_first_ionization = np.empty(0)
    std_distance_to_first_ionization = np.empty(0)
    
    frac_ioniz_in_Chodura = np.zeros(len(z_inner_target)-1)
    frac_ioniz_in_Debye = np.zeros(len(z_inner_target)-1)
    pindex_in_Chodura = np.empty(0,dtype=int)
    pindex_in_Debye = np.empty(0,dtype=int)
    frac_prompt_Chodura = np.zeros(len(z_inner_target)-1)
    frac_prompt_Debye = np.zeros(len(z_inner_target)-1)
    
    for seg in range(len(z_inner_target)-1):

        particle_index_list = np.empty(0,dtype=int)
        for p in range(int(nP)):
            
            if z[p,0]<=z_inner_target[seg] and z[p,0]>z_inner_target[seg+1]:
                particle_index_list = np.append(particle_index_list, int(p))
        
        distance_to_first_ionization = np.empty(0)
        Chodura_top, Debye_top = (0,0)
        for p in particle_index_list:
            
            time_index = np.nonzero(charge[p])[0]
            if time_index.size == 0:
                tally_never_ionizes[seg] += 1
            else:
                time_index = time_index[0]
                
                x1 = x[p,0]
                y1 = y[p,0]
                z1 = z[p,0]
                x2 = x[p,time_index]
                y2 = y[p,time_index]
                z2 = z[p,time_index]
                
                dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                distance_to_first_ionization = np.append(distance_to_first_ionization, dist)
                
                if dist<Chodura[seg]:
                    pindex_in_Chodura = np.append(pindex_in_Chodura, int(p))
                    Chodura_top += 1
                    
                    if dist<Debye[seg]:
                        pindex_in_Debye = np.append(pindex_in_Debye, int(p))
                        Debye_top += 1
                        
        avg_distance_to_first_ionization = np.append(avg_distance_to_first_ionization, \
                                                     np.average(distance_to_first_ionization))
        std_distance_to_first_ionization = np.append(std_distance_to_first_ionization, \
                                                     np.std(std_distance_to_first_ionization))
        frac_ioniz_in_Chodura[seg] = Chodura_top/len(particle_index_list)
        frac_ioniz_in_Debye[seg] = Debye_top/len(particle_index_list)
        
        if plotting[0]==1:
            plt.close()
            plt.hist(distance_to_first_ionization*1000,100) #in mm
            plt.xlabel('Distance [mm]')
            plt.ylabel('Counts')
            plt.show(block=True)
        
        ##############################################################################
        # Calculate fraction of particles ionizing in a sheath that promptly redeposit
        ##############################################################################
        
        angle_Chodura = angle[pindex_in_Chodura]
        angle_Debye = angle[pindex_in_Debye]
        
        is_prompt_Chodura = angle_Chodura <= 8.5
        is_prompt_Debye = angle_Debye <= 8.5
        
        frac_prompt_Chodura[seg] = np.sum(is_prompt_Chodura)/len(pindex_in_Chodura)
        frac_prompt_Debye[seg] = np.sum(is_prompt_Debye)/len(pindex_in_Debye)

    print('Never ionized:',tally_never_ionizes)
    
    if plotting[1]==1:
        plt.close()
        if tile_shift_indices != []:
            for i,v in enumerate(tile_shift_indices):
                if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
                else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
        if Bangle_shift_indices != []:
            for i,v in enumerate(Bangle_shift_indices):
                if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
                else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
        
        plt.errorbar(rmrsMid, avg_distance_to_first_ionization*1000, std_distance_to_first_ionization*1000)
        plt.scatter(rmrsMid, avg_distance_to_first_ionization*1000, s=15)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Distance [mm]')
        plt.title('Average Distance to First Ionization')
        plt.show(block=True)
        plt.savefig('plots/avg_dist_to_first_ioniz.png')
        
        plt.close()
        if tile_shift_indices != []:
            for i,v in enumerate(tile_shift_indices):
                if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
                else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
        if Bangle_shift_indices != []:
            for i,v in enumerate(Bangle_shift_indices):
                if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
                else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
        
        plt.plot(rmrsMid, frac_ioniz_in_Chodura, label='Chodura', color='orange')
        plt.plot(rmrsMid, frac_ioniz_in_Debye, label='Debye', color='red')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Fraction')
        plt.title('Fraction of Particles First Ionizing in the Sheath')
        plt.legend()
        plt.show(block=True)
        plt.savefig('plots/frac_ioniz_in_sheath.png')
        
        plt.close()
        if tile_shift_indices != []:
            for i,v in enumerate(tile_shift_indices):
                if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
                else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
        if Bangle_shift_indices != []:
            for i,v in enumerate(Bangle_shift_indices):
                if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\Psi_B$')
                else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
        
        plt.plot(rmrsMid, frac_prompt_Chodura, label='Chodura', color='orange')
        plt.plot(rmrsMid, frac_prompt_Debye, label='Debye', color='red')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Fraction')
        plt.title('Fraction of First Ionizations in Sheath that Redeposit')
        plt.legend()
        plt.show(block=True)
        plt.savefig('plots/frac_sheath_ioniz_prompt_redep.png')
    
    return

def prompt_redep_hist(inputs, fileDir, fileON, fileOFF):
    [nP10, dt10, nT10] = inputs
    
    pathON = fileDir+fileON
    pathOFF = fileDir+fileOFF
    
    posON = netCDF4.Dataset(pathON, "r", format="NETCDF4")
    posOFF = netCDF4.Dataset(pathOFF, "r", format="NETCDF4")
    
    angleON = posON.variables['angle'][:]/2/np.pi #radians -> rotations
    angleOFF = posOFF.variables['angle'][:]/2/np.pi #radians -> rotations
    timeON = posON.variables['time'][:]
    timeOFF = posOFF.variables['time'][:]
    
    is_promptON = angleON<=1
    is_promptOFF = angleOFF<=1
    
    frac_promptON = np.sum(is_promptON) / len(angleON)
    frac_promptOFF = np.sum(is_promptOFF) / len(angleOFF)
    
    print('Fraction promptly redeposited with SurfModel ON:', frac_promptON)
    print('Fraction promptly redeposited with SurfModel OFF:', frac_promptOFF)
    
    if 0:
        plt.close()
        bins = np.logspace(-3,4)
        #plt.hist(angleON, 200, (0,1), color='g', label='ON')
        #plt.hist(angleOFF, 200, (0,1), color='r', label='OFF')
        plt.axvline(x=1, color='k') #1 rotation = 2pi radians
        #setting density=True divides each bar by total counts and bin width to give a pdf
        plt.hist(angleON, bins, color='g', label='ON', alpha=0.5, density=False)
        plt.hist(angleOFF, bins, color='r', label='OFF', alpha=0.5, density=False)
        plt.xscale('log')
        plt.xlabel('Rotations')
        plt.ylabel('Counts')
        plt.title('Histogram of rotations completed before depositing \n with nP=1e'\
                  +str(nP10)+', dt=1e-'+str(dt10)+', nT=1e'+str(nT10))
        plt.legend()
    
    if 1: 
        plt.close()
        plt.hist(timeOFF, 200, (0,0.2), color='r', label='OFF')
        plt.hist(timeON, 200, (0,0.2), color='g', label='ON')
        plt.xlabel('Time [s]')
        plt.ylabel('Counts')
        plt.title('Flight time before striking \n nP=1e'+str(nP10))
        plt.legend()
    
    '''
    timeONzero = 10**nP10 - np.count_nonzero(timeON)
    timeOFFzero = 10**nP10 - np.count_nonzero(timeOFF)
    print('Particles that THINK they never leave the surface')
    print('Surface Model ON:', timeONzero)
    print('Surface Model OFF:', timeOFFzero)
    print('\n')
    '''

    return

if __name__ == "__main__":
    #plot_history2D('perlmutter/D3p5t8T5/history.nc')
    plot_surf_nc(100692963657457.89, 5, 5, [1,9], [3,8,9], "perlmutter/D3p5t8T5/surface.nc", norm=None)
    #analyze_leakage('perlmutter/history_D3t6.nc')
    #analyze_forces('gradT dv', 't', rzlim=True, colorbarLimits=[], dt=1e-8)
    
    #init()
    #plot_gitr_gridspace()
    #plot_particle_source()
    #plot_history2D('history-alpine.nc', plot_particle_source=1, markersize=2)
    #plot_history2D("../../../../GITR/scratch/output/history.nc")
    #plot_history2D("perlmutter/historyT4_dist_first_ioniz.nc")
    #plot_surf_nc(37914807680566.16, "/Users/Alyssa/Dev/SAS-VW-Data/netcdf_data/nP5/surf-5-6.nc", norm='C')
    #spectroscopy(3791480768056.615,specFile='specP6T6.nc')
    #ionization_analysis([0,1], 'perlmutter/historyT5_dist_first_ioniz.nc', 'perlmutter/positionsT5_dist_first_ioniz.nc', [1,9], [3,8,9], W_surf=np.arange(11,22))
    #prompt_redep_hist([5,9,4], 'perlmutter/p5t9T4/','positions_SurfModelON.nc','positions_SurfModelOFF.nc')

