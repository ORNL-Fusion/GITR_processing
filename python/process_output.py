import sys, os
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as path
import netCDF4
import solps

################################################
# setting directories and special constants
################################################

run_directory = '/Users/Alyssa/Dev/GITR/scratch'
#run_directory = '/pscratch/sd/h/hayes/sasvw-pa-fav-history'
setup_directory = '../examples/sasvw-pa-fav/setup'
rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt'

#prog angle
W_surf_indices = np.arange(11,22)
tile_shift_indices = [1,9]
Bangle_shift_indices = [3,8,9]
r_sp, z_sp = 1.49829829, 1.19672716 #prog angle & favorable
#r_sp, z_sp = 1.49829824, 1.19672712 #prog angle & unfavorable
'''
#vertex
W_surf_indices = np.arange(16,25)
tile_shift_indices = [2,6]
Bangle_shift_indices = [3,6]
r_sp, z_sp = 1.50230407, 1.23187366 #vertex & favorable
'''
sys.path.insert(0, os.path.abspath(setup_directory))
import makeParticleSource

################################################
# function definitions
################################################

def init(W_indices = W_surf_indices, plot_rz=False):
    profilesFile = setup_directory+'/../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    
    #check which target is the region of interest
    r_inner_target = profiles.variables['r_inner_target'][:]
    z_inner_target = profiles.variables['z_inner_target'][:]
    rmrs_coarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices]
    
    #plot r,z to check that target indices are restricted to W surfaces
    if plot_rz:
        plt.close()
        plt.plot(r_inner_target, z_inner_target)
    
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':2})
    plt.rcParams.update({'lines.markersize':1})

    return profiles, W_indices, r_inner_target, z_inner_target, rmrs_coarse

def init_geom(gitr_rz = setup_directory+'/assets/gitr_rz.txt', \
              rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt', \
              W_fine_file = setup_directory+'/assets/W_fine.txt'):
    
    #import refined rmrs at the W surface
    with open(rmrs_fine_file, 'r') as file:
        rmrs_fine = file.readlines()   
    rmrsFine = np.array(rmrs_fine,dtype='float')
    
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
    
    #R_surf = R[W_fine]
    #Z_surf = Z[W_fine]
    
    return W_fine, R, Z, rmrsFine

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

def plot_history2D(history_file, bFile=run_directory+'/input/bField.nc', \
                   basic=0, continuousChargeState=1, endChargeState=0, \
                   plot_particle_source=0, markersize=0):
    
    if plot_particle_source:
        particleSource = netCDF4.Dataset(run_directory+"/input/particleSource.nc", "r", format="NETCDF4")
        x0 = particleSource.variables['x'][:]
        z0 = particleSource.variables['z'][:]
    
    profiles, W_indices, R, Z, rmrs = init()
    W_fine, r_wall, z_wall, rmrs_fine = init_geom()
    r_target_fine = r_wall[W_fine]
    z_target_fine = z_wall[W_fine]
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")
    
    if bFile != None:
        bField = netCDF4.Dataset(bFile)
        r_bField = bField.variables['r'][:]
        z_bField = bField.variables['z'][:]
        psi = bField.variables['psi'][:]

    plt.rcParams.update({'lines.linewidth':0.3})
    plt.rcParams.update({'lines.markersize':markersize})
    plt.rcParams.update({'font.size':14})

    nP = len(history.dimensions['nP'])
    print('nP:',nP)
    nT = len(history.dimensions['nT'])
    x = history.variables['x'][:]
    y = history.variables['y'][:]
    z = history.variables['z'][:]
    r = np.sqrt(x**2 + y**2)
    charge = history.variables['charge'][:]

    plt.close()
    if plot_particle_source: plt.scatter(x0,z0,marker='o',s=10)
    plt.plot(r_wall, z_wall,'-k',linewidth=1.5)
    plt.plot(r_target_fine, z_target_fine,'-m',linewidth=2)
    plt.axis('scaled')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
        
    #define charge state to color mapping
    colors = {0:'black', 1:'firebrick', 2:'darkorange', 3:'gold', 4:'limegreen', 5:'dodgerblue', \
              6:'mediumpurple', 7:'darkviolet', 8:'darkmagenta', 9:'deeppink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            #if z[p][0]<=Z[W_indices][3] and z[p][0]>Z[W_indices][4]:
                print(p,'out of', nP)
                plt.plot(r[p][:],z[p][:])
                #plt.scatter(r[p][:],z[p][:],marker='o',s=5,c='b')
    
    counter=0
    if continuousChargeState==1:
        for p in np.arange(0,nP,1):
            t=0
            counter+=1
            while t<nT-1:
                if r[p][t] != r[p][t+1]: 
                    print("particle #", p, "moved at timestep", t, "with charge", charge[p][t])
                    plt.plot(r[p][t:t+2],z[p][t:t+2], colors[np.round(charge[p][t])])
                t+=1
    print('total particles:',counter)

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(r[p][:],z[p][:], colors[charge[p][-1]])
        plt.title('Particle Tracks by End Charge State')
    
    if bFile != None: plt.contour(r_bField,z_bField,psi,1,linestyles='--',colors='k',linewidths=1)
    #plt.scatter(1.49829829, 1.19672716, label='Strikepoint', marker='X', color='k', s=100, zorder=5)
    
    legend_dict = {'+0':'black', '+1':'firebrick', '+2':'darkorange', '+3':'gold', '+4':'limegreen', '+5':'dodgerblue', \
              '+6':'mediumpurple', '+7':'darkviolet', '+8':'darkmagenta', '+9':'deeppink', '+10':'gray'}
    
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    if basic==0: plt.legend(handles=patchList, fontsize=8, loc=2) #upper-left=2, lower-left=3
    
    #whole device
    #plt.xlim(1.0, 2.0)
    #plt.ylim(-1.5, 1.5)
    #pa-fav
    plt.xlim(1.37, 1.52)
    plt.ylim(1.06, 1.23)
    #vertex-fav
    #plt.xlim(1.0, 1.53)
    #plt.ylim(1.0, 1.23)
    #vertex-unfav
    #plt.xlim(1.0, 1.53)
    #plt.ylim(1.0, 1.23)
    
    plt.title('W Trajectories', fontsize=24)
    #plt.show(block=False)
    plt.savefig(run_directory+'/output/history.svg')
    plt.close()
    return

def plot_surf_nc(nP10, dt10, nT10, \
                 surface_file="surface.nc", positions_file='positions.nc', \
                 gitr_rz=setup_directory+'/assets/gitr_rz.txt', \
                 W_fine_file=setup_directory+'/assets/W_fine.txt', \
                 rmrs_fine_file=setup_directory+'/assets/rmrs_fine.txt', \
                 norm=None, plot_cumsum=0):
    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init()
    rmrsCoords = profiles.variables['rmrs_inner_target'][W_indices]
    surface = netCDF4.Dataset(surface_file, "r", format="NETCDF4")
    #pps_per_nP = 1394457587.59747
    
    pps_per_nP, partSource_flux, fluxD, fluxC = makeParticleSource.distributed_source(nP=(nP10[0] * (10**int(nP10[1]))), \
                surfW = W_surf_indices, \
                tile_shift_indices = tile_shift_indices, \
                Bangle_shift_indices = Bangle_shift_indices, \
                geom = setup_directory+'/../input/gitrGeometry.cfg', \
                profiles_file = setup_directory+'/../input/plasmaProfiles.nc', \
                gitr_rz = setup_directory+'/assets/gitr_rz.txt', \
                rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt', \
                W_fine_file = setup_directory+'/assets/W_fine.txt', \
                ftDFile = setup_directory+'/assets/ftridynBackgroundD.nc', \
                ftCFile = setup_directory+'/assets/ftridynBackgroundC.nc', \
                ftWFile = setup_directory+'/../input/ftridynSelf.nc', \
                configuration = 'random', \
                use_fractal_tridyn_outgoing_IEADS = 1, \
                plot_variables = 0)
    
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
    
    if positions_file != '':
        positions = netCDF4.Dataset(positions_file, "r", format="NETCDF4")
        flightAngle = positions.variables['angle'][:]
        is_prompt_redep = flightAngle < 2*np.pi
        prompt_redep_rate = np.sum(is_prompt_redep) / len(flightAngle)

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
    print('total gross eroded flux',sum(grossEro*area)/sum(area))
    print('total redeposited flux',sum(grossDep*area)/sum(area))
    print('total net deposited flux',sum(netDep*area)/sum(area))
    print('redeposition rate',100 * (sum(grossDep*area)/sum(area)) / (sum(grossEro*area)/sum(area)), '%')
    print('prompt redeposition rate', 100 * prompt_redep_rate, '%')
    #print('self-sputtering fraction',100 * (sum(grossEro)-sum(partSource_flux))/sum(grossEro), '%')
    '''
    if norm=='C':
        grossEro_norm = grossEro/fluxC
        grossDep_norm = grossDep/fluxC
        netDep_norm = netDep/fluxC
    elif norm=='D':
        grossEro_norm = grossEro/fluxD
        grossDep_norm = grossDep/fluxD
        netDep_norm = netDep/fluxD
    else: '''
    grossEro_norm = grossEro
    grossDep_norm = grossDep
    netDep_norm = netDep
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netDep_cumsum = np.cumsum(netDep)
    
    #take gross erosion of a slice in the filterscope range
    r1_start, z1_start = 1.491288, 1.21629
    r1_end, z1_end = 1.49274, 1.22149
    #account for when View 1 is before or after OSP
    is_neg = -1
    if 'vertex' in setup_directory:
        r1_start, z1_start = 1.49274, 1.22149
        r1_end, z1_end = 1.491288, 1.21629
        is_neg=1
    rmrs1_start = is_neg*np.sqrt((z1_start-z_sp)**2 + (r1_start-r_sp)**2)
    rmrs1_end = is_neg*np.sqrt((z1_end-z_sp)**2 + (r1_end-r_sp)**2)
    
    r2_start, z2_start = 1.492041, 1.19418 
    r2_end, z2_end = 1.492716, 1.18709 
    rmrs2_start = np.sqrt((z2_start-z_sp)**2 + (r2_start-r_sp)**2)
    rmrs2_end = np.sqrt((z2_end-z_sp)**2 + (r2_end-r_sp)**2)
    
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
    
    if 'vertex' in setup_directory:
        lower_bound = rmrs1_start
        upper_bound = rmrs1_end
    else:
        lower_bound = rmrs1_end
        upper_bound = rmrs1_start
        
    for i,v in enumerate(rmrsFine):
        if v >= lower_bound and v <= upper_bound:
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
    if 'vertex' in setup_directory:
        V1_grossEro += grossEro[V1_indices[0]-1] * (rmrs1_start-rmrsFine[V1_indices[0]])
    else:
        V1_grossEro += grossEro[V1_indices[0]-1] * (rmrs1_end-rmrsFine[V1_indices[0]])
    V1_area += rmrs1_start-rmrsFine[V1_indices[0]]
    V2_grossEro += grossEro[V2_indices[0]-1] * (rmrs2_start-rmrsFine[V2_indices[0]])
    V2_area += rmrs2_start-rmrsFine[V2_indices[0]]
    V3_grossEro += grossEro[V3_indices[0]-1] * (rmrs3_start-rmrsFine[V3_indices[0]])
    V3_area += rmrs3_start-rmrsFine[V3_indices[0]]

    #subtract fraction of last cell starting inside the view
    if 'vertex' in setup_directory:
        V1_grossEro -= grossEro[V1_indices[-1]] * (rmrs1_end-rmrsFine[V1_indices[-1]+1])
    else:
        V1_grossEro -= grossEro[V1_indices[-1]] * (rmrs1_start-rmrsFine[V1_indices[-1]+1])
    V1_area -= rmrs1_start-rmrsFine[V1_indices[-1]+1] 
    V2_grossEro -= grossEro[V2_indices[-1]] * (rmrs2_end-rmrsFine[V2_indices[-1]+1])
    V2_area -= rmrs2_end-rmrsFine[V2_indices[-1]+1]
    V3_grossEro -= grossEro[V3_indices[-1]] * (rmrs3_end-rmrsFine[V3_indices[-1]+1])
    V3_area -= rmrs3_end-rmrsFine[V3_indices[-1]+1]
    
    V1_grossEro = V1_grossEro / V1_area
    V2_grossEro = V2_grossEro / V2_area
    V3_grossEro = V3_grossEro / V3_area
    
    #filterscope_diameter1 = 0.005404199 #option to replace V1_area in the toroidal_fraction1
    tokamak_circumference1 = 2 * np.pi * 1.492014
    toroidal_fraction1 = V1_area / tokamak_circumference1
    view_fraction1 = toroidal_fraction1 * np.pi/4

    #filterscope_diameter2 = 0.007125431
    tokamak_circumference2 = 2 * np.pi * 1.4923785
    toroidal_fraction2 = V2_area / tokamak_circumference2
    view_fraction2 = toroidal_fraction2 * np.pi/4

    #filterscope_diameter3 = 0.007476819
    tokamak_circumference3 = 2 * np.pi * 1.492046
    toroidal_fraction3 = V3_area / tokamak_circumference3
    view_fraction3 = toroidal_fraction3 * np.pi/4
    
    V1_grossEro = view_fraction1 * V1_grossEro
    V2_grossEro = view_fraction2 * V2_grossEro
    V3_grossEro = view_fraction3 * V3_grossEro
    
    print('\n')
    print('gross erosion in View 1:', V1_grossEro)
    print('gross erosion in View 2:', V2_grossEro)
    print('gross erosion in View 3:', V3_grossEro)
    
    plt.rcParams.update({'font.size':14})
    plt.rcParams.update({'lines.linewidth':5}) 
    '''
    #plot self-sputtering
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsFine, 100*(grossEro-partSource_flux)/grossEro)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Percentage')
    plt.title('Percentage of Gross Erosion from Self-Sputtering')
    plt.show(block=True)
    '''
    #plot main surface plot with 3 views
    plt.close()
    plt.axvspan(rmrs1_start, rmrs1_end, color='#f7bc00', alpha=0.5)
    plt.axvspan(rmrs2_start, rmrs2_end, color='lightsalmon', alpha=0.5)
    plt.axvspan(rmrs3_start, rmrs3_end, color='#f99301', alpha=0.5)
    
    plt.plot(rmrsFine,np.zeros(len(rmrsFine)),'gray')
    '''
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    '''
    plt.plot(rmrsFine,grossEro_norm,'r', label='Gross Erosion')
    plt.plot(rmrsFine,grossDep_norm,'g', label='Gross Deposition')
    plt.plot(rmrsFine,netDep_norm,'k', label='Net Deposition')
    
    plt.xlabel('D-Dsep [m]')
    if norm!=None: plt.ylabel('\u0393$_{W,outgoing}$ / \u0393$_{%s,incoming}$'%norm)
    if norm==None: plt.ylabel('\u0393$_{W}$ [m$^{-2}$s$^{-1}$]')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    plt.legend(fontsize=10)#loc='upper left')
    plt.title('GITR-Predicted Erosion and Deposition', fontsize=20)
    #plt.title('GITR Predicted Erosion and Redeposition Profiles,\nnP='+str(nP10[0])+'e'+str(nP10[1])+', dt=1e-'+str(dt10)+', nT='+str(nT10[0])+'e'+str(nT10[1]), fontsize=30)
    plt.show(block=True)
    #plt.savefig('plots/surface.png')
    
    if plot_cumsum:
        plt.close()
        plt.plot(rmrsFine,grossEro_cumsum,'r', label='Gross Erosion')
        plt.plot(rmrsFine,grossDep_cumsum,'g', label='Gross Deposition')
        plt.plot(rmrsFine,netDep_cumsum,'k', label='Net Erosion')
        plt.yscale('log')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Flux [m$^{-2}$s$^{-1}$]')
        plt.legend(loc='upper left')
        plt.title('GITR Predicted Cumulative Sum of \n Erosion and Redeposition Profiles',fontsize=20)
        plt.savefig('plots/surface_cumsum.pdf')
    
    return
        
def spectroscopy(pps_per_nP, View=3, \
                 specFile='spec.nc',plotting=1):
    
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
    #density_neutrals_unitless = density_unitless[0]
    #temporarily make this a summation over all W charge states
    density_neutrals_unitless = np.sum(density_unitless,axis=0) 
    
    plt.pcolor(gridr,gridz,density_neutrals_unitless)
    plt.axis('Scaled')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.colorbar(label='# Computational Neutrals')
    plt.title('Raw W0 Spec Output form GITR')
    plt.savefig('plots/spec_compNeutrals.png')
    
    if View==1:
        filterscope_diameter = 0.005404199
        tokamak_circumference = 2 * np.pi * 1.492014
    elif View==2:
        filterscope_diameter = 0.007125431
        tokamak_circumference = 2 * np.pi * 1.4923785
    elif View==3:
        filterscope_diameter = 0.007476819
        tokamak_circumference = 2 * np.pi * 1.492046
    
    toroidal_fraction = filterscope_diameter / tokamak_circumference
    view_fraction = toroidal_fraction * np.pi/4
    
    dA = np.pi * (filterscope_diameter/2)**2
    neutral_flux = density_neutrals_unitless * pps_per_nP / dA # W0 m-2 s-1
    neutral_flux_sliced = view_fraction * neutral_flux

    '''
    r_mid = dr/2 + rr[:-1,:]
    dV = 2 * np.pi * r_mid[:,:-1] * dr * dz * toroidal_fraction
    density_volumetric = pps_per_nP * (density_neutrals/dV) * 1e-8 * 3e3 # W0 m-2 s-1
    #neutral_flux = density_volumetric * 3e3 # W0 m-2 s-1

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
    '''
    if View==1:
        line1 = -0.564378*gridr + 2.063959
        line2 = -0.790843*gridr + 2.395664
        
        if plotting:
            plt.close()
            plt.plot(R,Z)
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,neutral_flux_sliced)
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('r [m]')
            plt.ylabel('z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
            plt.title('Toroidal Slice of W0 Density')
            plt.savefig('plots/spec_density_sliced.png')
        
        rstart, rend = 62, 71
        zstart, zend = 95, 104
        
        gridr = gridr[rstart:rend]
        gridz = gridz[zstart:zend]
        fscope = neutral_flux_sliced[zstart:zend, rstart:rend]
        
        line1 = -0.564378*gridr + 2.063959
        line2 = -0.790843*gridr + 2.395664
        
        fscope[np.where(fscope==0)] = 'NaN'
        fscope[-1,1] = fscope[-1,5:] = fscope[-2,6:] = fscope[-3:,7:] = 'NaN'
        for i in range(0,8):
            fscope[i,:-i-2] = 'NaN'
        
    elif View==2: 
        line1 = -0.685413*gridr + 2.216844
        line2 = -0.744624*gridr + 2.298601
        
        if plotting:
            plt.close()
            plt.plot(R,Z)
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,neutral_flux_sliced)
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('R [m]')
            plt.ylabel('Z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
            plt.title('Toroidal Slice of W0 Density')
            plt.savefig('plots/spec_density_sliced.png')
        
        rstart, rend = 50, 72
        zstart, zend = 69, 92
        
        gridr = gridr[rstart:rend]
        gridz = gridz[zstart:zend]
        fscope = neutral_flux_sliced[zstart:zend, rstart:rend]
        
        line1 = -0.685413*gridr + 2.216844
        line2 = -0.744624*gridr + 2.298601
        
        fscope[np.where(fscope==0)] = 'NaN'
        fscope[0] = fscope[-1] = fscope[:,0] = 'NaN'
        for i in range(1,7):
            fscope[-i,i+2:] = 'NaN'
        for i in range(7,9):
            fscope[-i,i+3:] = 'NaN'
        for i in range(9,18):
            fscope[-i,i+4:] = 'NaN'
        for i in range(7,12):
            fscope[-i,i+3:] = 'NaN'
        for i in range(0,11):
            fscope[:-i-3,i] = 'NaN'
        for i in range(11,20):
            fscope[:-i-2,i] = 'NaN'
        
    elif View==3:     
        line1 = -0.906176*gridr + 2.515464
        line2 = -0.972056*gridr + 2.604085
        
        if plotting:
            plt.close()
            plt.plot(R,Z)
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,neutral_flux_sliced)
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('R [m]')
            plt.ylabel('Z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
            plt.title('Toroidal Slice of W0 Density')
            plt.savefig('plots/spec_density_sliced.png')
        
        rstart, rend = 39, 72
        zstart, zend = 40, 80
        
        gridr = gridr[rstart:rend]
        gridz = gridz[zstart:zend]
        fscope = neutral_flux_sliced[zstart:zend, rstart:rend]
        
        line1 = -0.906176*gridr + 2.515464
        line2 = -0.972056*gridr + 2.604085
        
        fscope[np.where(fscope==0)] = 'NaN'
        for i in range(0,35):
            fscope[-i,i+2:] = 'NaN'
        for i in range(0,40):
            fscope[i,:-1-i] = 'NaN'
        for i in range(0,20):
            fscope[-i-4,:i] = 'NaN'
        for i in range(0,16):
            fscope[24-i,17+i] = 'NaN'
        fscope[9:,-2:] = fscope[11:,-4:] = 'NaN'
        fscope[12,28] = fscope[13,19] = fscope[12,20] = 'NaN'
    
    if True:
        plt.close()
        plt.plot(R,Z)
        plt.plot(gridr,line1,'orange')
        plt.plot(gridr,line2,'orange')
        plt.pcolor(gridr,gridz,fscope,shading='auto')
        plt.axis('Scaled')
        plt.xlim(gridr[0],gridr[-1])
        plt.ylim(gridz[0],gridz[-1])
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.colorbar(label='\n Gross Eroded W Flux [m$^{-2}$ s$^{-1}$]') #'\n Density [m$^{-3}$ s$^{-1}$]')

        fscope[np.isnan(fscope)] = 0
        print('W0 Neutral Flux:', np.sum(fscope)) # in units of m-2 s-1
        plt.title('Toroidal Slice of W0 Density \n View '+str(View)+' W0: %10e m$^{-2}$s$^{-1}$' %np.sum(fscope))
        plt.savefig('plots/spec_filterscope.png')
    
def analyze_leakage(historyFile, bFile = '../input/bField.nc'):
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
    gitr_rz=setup_directory+'/assets/gitr_rz.txt'
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
        
    r_wall = np.zeros(len(wall))
    z_wall = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        r_wall[i] = float(point[0])
        z_wall[i] = float(point[1])
    
    profilesFile = setup_directory+'/../input/plasmaProfiles.nc'
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
    plt.rcParams.update({'image.cmap':'coolwarm'})
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
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
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

def ionization_analysis(plotting, output_dir, historyFile, positionsFile, \
                        use_coarse_surfs=0, \
                        gitr_rz=setup_directory+'/assets/gitr_rz.txt', \
                        W_fine_file=setup_directory+'/assets/W_fine.txt', \
                        rmrs_fine_file=setup_directory+'/assets/rmrs_fine.txt'):    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init(W_surf_indices)
    history = netCDF4.Dataset(output_dir+historyFile)
    positions = netCDF4.Dataset(output_dir+positionsFile)
    
    z_inner_target = z_inner_target[W_surf_indices]
    rmrsMid = rmrs[:-1]
    rmrsCoords = profiles.variables['rmrs_inner_target'][W_surf_indices]
    nP = len(history.dimensions['nP'])
    charge = history.variables['charge'][:]
    perpDistanceToSurface = history.variables['perpDistanceToSurface'][:]
    angle = positions.variables['angle'][:] #final amount of radians travelled at the end of a particle packet's life
    hitWall = positions.variables['hitWall'][:]
    
    Debye, Chodura = theoretical_sheath(profiles,W_surf_indices)
    
    if use_coarse_surfs:
        #TODO: this still uses distance travelled and NOT perpendicular distance to surface
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
        
        rmrsPlotting = rmrsMid
        
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

    else:
        z0 = history.variables['z'][:][:,0] #just take the first timestep of each particle
        
        #import wall geometry to plot over fine surface
        with open(gitr_rz, 'r') as file:
            wall = file.readlines()
        
        #import W surface indices
        with open(W_fine_file, 'r') as file:
            W_fine = file.readlines()
        W_fine = np.array(W_fine,dtype='int')
        
        Z = np.zeros(len(wall))
        for i,line in enumerate(wall):
            point = line.split()
            Z[i] = float(point[1])
        
        #import refined rmrs at the W surface
        with open(rmrs_fine_file, 'r') as file:
            rmrs_fine = file.readlines()   
        rmrsFine = np.array(rmrs_fine,dtype='float')
        rmrsPlotting = rmrsFine
        
        Z = Z[W_fine]
        
        Z1 = Z[:-1]
        Z2 = Z[1:]
        
        tally_never_ionizes = np.zeros(len(Z1))
        particles_per_seg = np.zeros(len(Z1))
        avg_distance_to_first_ionization = np.zeros(len(Z1))
        std_distance_to_first_ionization = np.zeros(len(Z1))
        median_distance_to_first_ionization = np.zeros(len(Z1))
        
        frac_ioniz_in_Chodura = np.zeros(len(Z1))
        frac_ioniz_in_Debye = np.zeros(len(Z1))
        frac_prompt_Chodura = np.zeros(len(Z1))
        frac_prompt_Debye = np.zeros(len(Z1))
        frac_noHitWall_Chodura = np.zeros(len(Z1))
        frac_noHitWall_Debye = np.zeros(len(Z1))
        
        z_coarse_index = 0
        print('COARSE SEGMENT: 0')
        for seg in range(len(Z1)):
            
            if Z2[seg] < z_inner_target[z_coarse_index+1]:
                z_coarse_index += 1
                print('\n COARSE SEGMENT:', z_coarse_index)
    
            particle_index_list = np.empty(0,dtype=int)
            for p in range(int(nP)):
                
                if z0[p]<Z1[seg] and z0[p]>=Z2[seg]:
                    particle_index_list = np.append(particle_index_list, int(p))
            particles_per_seg[seg] = len(particle_index_list)
            
            distance_to_first_ionization = np.empty(0)
            pindex_in_Chodura = np.empty(0,dtype=int)
            pindex_in_Debye = np.empty(0,dtype=int)
            Chodura_top, Debye_top = (0,0)
            for p in particle_index_list:
                
                time_index = np.nonzero(charge[p])[0]
                if time_index.size == 0:
                    tally_never_ionizes[seg] += 1
                else:
                    time_index = time_index[0]
                    '''
                    x1 = x[p,0]
                    y1 = y[p,0]
                    z1 = z[p,0]
                    x2 = x[p,time_index]
                    y2 = y[p,time_index]
                    z2 = z[p,time_index]
                    
                    dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)'''
                    dist = perpDistanceToSurface[p,time_index]
                    distance_to_first_ionization = np.append(distance_to_first_ionization, dist)
                    
                    if dist<Chodura[z_coarse_index]:
                        pindex_in_Chodura = np.append(pindex_in_Chodura, int(p))
                        Chodura_top += 1
                        
                        if dist<Debye[z_coarse_index]:
                            pindex_in_Debye = np.append(pindex_in_Debye, int(p))
                            Debye_top += 1
            
            if len(distance_to_first_ionization) == 0:
                avg_distance_to_first_ionization[seg] = 0
                std_distance_to_first_ionization[seg] = 0
                median_distance_to_first_ionization[seg] = 0
            else:
                avg_distance_to_first_ionization[seg] = np.average(distance_to_first_ionization)
                std_distance_to_first_ionization[seg] = np.std(distance_to_first_ionization)
                median_distance_to_first_ionization[seg] = np.median(distance_to_first_ionization)
            
            if len(particle_index_list) == 0:
                frac_ioniz_in_Chodura[seg] = 0
                frac_ioniz_in_Debye[seg] = 0
            else:
                frac_ioniz_in_Chodura[seg] = Chodura_top/len(particle_index_list)
                frac_ioniz_in_Debye[seg] = Debye_top/len(particle_index_list)
            
            ##############################################################################
            # Calculate fraction of particles ionizing in a sheath that promptly redeposit
            ##############################################################################
            
            if len(pindex_in_Chodura) == 0:
                frac_prompt_Chodura[seg] = 0
                frac_noHitWall_Chodura[seg] = 0
            else:
                angle_Chodura = angle[pindex_in_Chodura]                
                is_prompt_Chodura = angle_Chodura <= 2*np.pi                
                frac_prompt_Chodura[seg] = np.sum(is_prompt_Chodura)/len(pindex_in_Chodura)
                
                hitWall_Chodura = hitWall[pindex_in_Chodura]
                frac_hitWall_Chodura = np.sum(hitWall_Chodura)/len(pindex_in_Chodura)
                frac_noHitWall_Chodura[seg] = 1-frac_hitWall_Chodura
    
            if len(pindex_in_Debye) == 0:
                frac_prompt_Debye[seg] = 0
                frac_noHitWall_Debye[seg] = 0
            else:
                angle_Debye = angle[pindex_in_Debye]
                is_prompt_Debye = angle_Debye <= 2*np.pi
                frac_prompt_Debye[seg] = np.sum(is_prompt_Debye)/len(pindex_in_Debye)
                
                hitWall_Debye = hitWall[pindex_in_Debye]
                frac_hitWall_Debye = np.sum(hitWall_Debye)/len(pindex_in_Debye)
                frac_noHitWall_Debye[seg] = 1-frac_hitWall_Debye
    
            #print('Segment #' + str(seg) + ': ' + tally_never_ionizes ' never ionizes of ' + str(particles_per_seg) ' total particles')
            print('Segment:', seg)
            print('Never ionized:',tally_never_ionizes[seg])
            print('Particles so far:', np.sum(particles_per_seg[:seg+1]))
            
            if plotting[0]==1:
                plt.close()
                plt.hist(distance_to_first_ionization*1000,100) #in mm
                plt.xlabel('Distance [mm]')
                plt.ylabel('Counts')
                plt.title('Histogram of Distance to 1st Ionization \n Segment: '+str(seg))
                plt.show(block=True)
    
    #############################################################################
    # Plotting ### Plotting ### Plotting ### Plotting ### Plotting ### Plotting #
    #############################################################################
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.errorbar(rmrsPlotting, avg_distance_to_first_ionization*1000, std_distance_to_first_ionization*1000, ecolor='lightpink', fmt='saddlebrown')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Distance [mm]')
    plt.title('Average Distance to First Ionization')
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_avg.png')
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
            
    plt.plot(rmrsPlotting, median_distance_to_first_ionization*1000, marker='o', color='saddlebrown')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Distance [mm]')
    plt.title('Median Distance to First Ionization')
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_med.png')
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting[31:], frac_ioniz_in_Chodura[31:], label='Chodura', color='darkorange')
    plt.plot(rmrsPlotting, frac_ioniz_in_Debye, label='Debye', color='crimson')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Fraction')
    plt.title('Fraction of Particles First Ionizing in the Sheath')
    plt.legend()
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_sheath.png')
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting, frac_prompt_Chodura, label='Chodura', color='darkorange')
    plt.plot(rmrsPlotting, frac_prompt_Debye, label='Debye', color='crimson')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Fraction')
    plt.title('Fraction of First Ionizations in Sheath \n that Promptly Redeposit')
    plt.legend()
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_prompt.png')
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting, frac_noHitWall_Chodura, label='Chodura', color='darkorange')
    plt.plot(rmrsPlotting, frac_noHitWall_Debye, label='Debye', color='crimson')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Fraction')
    plt.title('Fraction of First Ionizations in Sheath \n that Never Hit a Wall')
    plt.legend()
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_nohitwall.png')
    
    frac_delayed_Chodura = 1 - frac_prompt_Chodura - frac_noHitWall_Chodura
    frac_delayed_Debye = 1 - frac_prompt_Debye - frac_noHitWall_Debye
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
            
    plt.plot(rmrsPlotting, frac_prompt_Chodura, label='Prompt Redeposition', color='saddlebrown')
    plt.plot(rmrsPlotting, frac_delayed_Chodura, label='Delayed Redeposition', color='darkorange')
    plt.plot(rmrsPlotting, frac_noHitWall_Chodura, label='Never Hit a Wall', color='burlywood')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Fraction')
    plt.title('Behavior of Particles that First Ionize \n in the Chodura Sheath')
    plt.legend()
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_Chodura.png')
    rmrsMidpoints = rmrsPlotting[1:] + (rmrsPlotting[1:]-rmrsPlotting[:-1])/2
    top = frac_prompt_Chodura[1:-1] * (rmrsMidpoints[1:]-rmrsMidpoints[:-1])
    bottom = rmrsMidpoints[-1]-rmrsMidpoints[1]
    print('Unweighted Average Prompt Fraction (Chodura):', np.average(frac_prompt_Chodura))
    print('Weighted Average Prompt Fraction (Chodura):', np.sum(top) / bottom)
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting, frac_prompt_Debye, label='Prompt Redeposition', color='maroon')
    plt.plot(rmrsPlotting, frac_delayed_Debye, label='Delayed Redeposition', color='crimson')
    plt.plot(rmrsPlotting, frac_noHitWall_Debye, label='Never Hit a Wall', color='lightcoral')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Fraction')
    plt.title('Behavior of Particles that First Ionize \n in the Debye Sheath')
    plt.legend()
    if plotting[1]==1: plt.show(block=True)
    else: plt.savefig(output_dir+'plots/ioniz_Debye.png')
    
    overall_avg = np.average(avg_distance_to_first_ionization[:])
    leg2_avg = np.average(avg_distance_to_first_ionization[15:70])
    leg3_avg = np.average(avg_distance_to_first_ionization[70:])
    print('\nOverall average perpendicular distance to first ionization:\n', overall_avg, 'meters')
    print('\nLeg 2 average perpendicular distance to first ionization:\n', leg2_avg, 'meters')
    print('\nLeg3 average perpendicular distance to first ionization:\n', leg3_avg, 'meters')
    
    print('\nFraction that first ionize in the Debye sheath past OSP:')
    print('Min:', np.min(frac_ioniz_in_Debye[32:]))
    print('Max:', np.max(frac_ioniz_in_Debye[31:]))
    print('Average:', np.average(frac_ioniz_in_Debye[31:]))
    
    print('\nFraction that first ionize in the Chodura sheath past OSP:')
    print('Min:', np.min(frac_ioniz_in_Chodura[31:]))
    print('Max:', np.max(frac_ioniz_in_Chodura[31:]))
    print('Average:', np.average(frac_ioniz_in_Chodura[31:]))  
    
    print('\nFraction of first ionizations in the Chodura sheath that promptly redeposit past OSP:')
    print('Min:', np.min(frac_prompt_Chodura[31:]))
    print('Max:', np.max(frac_prompt_Chodura[31:]))
    print('Average:', np.average(frac_prompt_Chodura[31:])) 
    return

def prompt_redep_hist(inputs, fileDir, fileOFF, fileON=None):
    [nP10, dt10, nT10] = inputs
    
    pathOFF = fileDir+fileOFF
    posOFF = netCDF4.Dataset(pathOFF, "r", format="NETCDF4")
    angleOFF = posOFF.variables['angle'][:]/2/np.pi #radians -> rotations
    timeOFF = posOFF.variables['time'][:]
    is_promptOFF = angleOFF<=1
    frac_promptOFF = np.sum(is_promptOFF) / len(angleOFF)
    print('Fraction promptly redeposited with SurfModel OFF:', frac_promptOFF)
    
    if fileON != None: 
        pathON = fileDir+fileON
        posON = netCDF4.Dataset(pathON, "r", format="NETCDF4")
        angleON = posON.variables['angle'][:]/2/np.pi #radians -> rotations
        timeON = posON.variables['time'][:]
        is_promptON = angleON<=1
        frac_promptON = np.sum(is_promptON) / len(angleON)
        print('Fraction promptly redeposited with SurfModel ON:', frac_promptON)
    
    if 1:
        plt.close()
        bins = np.logspace(-3,4)
        #plt.hist(angleON, 200, (0,1), color='g', label='ON')
        #plt.hist(angleOFF, 200, (0,1), color='r', label='OFF')
        plt.axvline(x=1, color='k') #1 rotation = 2pi radians
        #setting density=True divides each bar by total counts and bin width to give a pdf
        if fileON != None: plt.hist(angleON, bins, color='g', label='ON', alpha=0.5, density=False)
        plt.hist(angleOFF, bins, color='r', label='OFF', alpha=0.5, density=False)
        plt.xscale('log')
        plt.xlabel('Rotations')
        plt.ylabel('Counts')
        plt.title('Histogram of rotations completed before depositing \n with nP=1e'\
                  +str(nP10)+', dt=1e-'+str(dt10)+', nT=1e'+str(nT10))
        plt.legend()
    
    if 0: 
        plt.close()
        plt.hist(timeOFF, 200, (0,0.2), color='r', label='OFF')
        if fileON != None: plt.hist(timeON, 200, (0,0.2), color='g', label='ON')
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

def particle_diagnostics_hist(nP_input, pdFile, segment_counter=50, seg_hist_plotting=0, plot_blocker=1):
    
    diagnostics = netCDF4.Dataset(pdFile, "r", format="NETCDF4")
    bin_edges_time = diagnostics.variables['bin_edges_time'][:]
    histogram_particle_time = diagnostics.variables['histogram_particle_time'][:][:-1]
    bin_width_time = bin_edges_time[1]-bin_edges_time[0]
    bin_edges_angle = diagnostics.variables['bin_edges_angle'][:]
    histogram_particle_angle = diagnostics.variables['histogram_particle_angle'][:][:-1]
    bin_width_angle = bin_edges_angle[1]-bin_edges_angle[0]
    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init()
    rmrsCoords = profiles.variables['rmrs_inner_target'][W_surf_indices]
    
    #import refined rmrs at the W surface
    with open(rmrs_fine_file, 'r') as file:
        rmrs_fine = file.readlines()   
    rmrsFine = np.array(rmrs_fine,dtype='float')
    rmrsPlotting = rmrsFine
    
    ####################################################################################
    # Option for plotting flight time and flight angle histograms over 1 surface segment
    ####################################################################################
    
    if seg_hist_plotting:
        histogram_particle_time_real = histogram_particle_time
        histogram_particle_angle_real = histogram_particle_angle
        nP = np.sum(histogram_particle_time_real[segment_counter])
        
        plt.close()
        plt.bar(bin_edges_time[:-1], histogram_particle_time_real[segment_counter], width=bin_width_time, color='darkviolet', align='edge', edgecolor='k')
        plt.xlabel('Logarithmic Time [log(sec)]')
        plt.ylabel('Counts\n')
        plt.title('Flight Time before Striking the Surface \n Surface %i, nP=%.4E'%(segment_counter,nP))
        plt.show(block=False)
        
        plt.close()
        plt.bar(bin_edges_angle[:-1], histogram_particle_angle_real[segment_counter], width=bin_width_angle, color='darkkhaki', align='edge', edgecolor='k')
        plt.xlabel('Angle [rad]')
        plt.ylabel('Counts')
        plt.title('Flight Angle before Striking the Surface \n Surface %i, nP=%.4E'%(segment_counter,nP))
        plt.show(block=False)
        
    ###################################################################
    # Calculate average, std, and median flight times and flight angles
    ###################################################################
    
    num_segments = len(histogram_particle_time)
    avg_flight_time = np.zeros(num_segments)
    std_flight_time = np.zeros(num_segments)
    med_flight_time = np.zeros(num_segments)
    avg_flight_angle = np.zeros(num_segments)
    std_flight_angle = np.zeros(num_segments)
    med_flight_angle = np.zeros(num_segments)
    
    for seg in range(num_segments):
        if np.sum(histogram_particle_time[seg]) == 0:
            avg_flight_time[seg] = 0
            std_flight_time[seg] = 0
            med_flight_time[seg] = 0
        else:
            avg_flight_time[seg] = np.sum((histogram_particle_time[seg]*bin_edges_time[:-1]))/np.sum(histogram_particle_time[seg])
            
            time_bin_has_particles_binned = histogram_particle_time[seg]>0 #bool to zero out time bins with no particles
            med_flight_time[seg] = np.median(time_bin_has_particles_binned * bin_edges_time[:-1])
        
        if np.sum(histogram_particle_angle[seg]) == 0:
            avg_flight_angle[seg] = 0
            std_flight_angle[seg] = 0
            med_flight_angle[seg] = 0
        else:
            avg_flight_angle[seg] = np.sum((histogram_particle_angle[seg]*bin_edges_angle[:-1]))/np.sum(histogram_particle_angle[seg])
            
            angle_bin_has_particles_binned = histogram_particle_angle[seg]>0 #bool to zero out angle bins with no particles
            med_flight_angle[seg] = np.median(angle_bin_has_particles_binned * bin_edges_angle[:-1])
    
    avg_flight_time[np.where(avg_flight_time==0)] = 'NaN'
    avg_flight_time = 10**avg_flight_time
    med_flight_time = 10**med_flight_time
    #print('\n'+'Median Flight Time before Striking the Surface:', med_flight_time, '\n')
    
    avg_flight_angle[np.where(avg_flight_angle==0)] = 'NaN'
    avg_flight_angle = 10**avg_flight_angle
    med_flight_angle = 10**med_flight_angle
    #print('\n'+'Median Flight Angle before Striking the Surface:', med_flight_angle, '\n')
    
    ###########
    # Plotting
    ###########
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting, avg_flight_time, 'darkviolet')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flight Time [s]')
    plt.title('Average Flight Time before Striking the Surface')
    #plt.legend()
    plt.show(block=plot_blocker)
    #if not plot_blocker: plt.savefig('plots/avg_flight_time.png')
    
    plt.close()
    if tile_shift_indices != []:
        for i,v in enumerate(tile_shift_indices):
            if i==0: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed', label='Walll\nVertices')
            else: plt.axvline(x=rmrsCoords[v], color='k', linestyle='dashed')
    if Bangle_shift_indices != []:
        for i,v in enumerate(Bangle_shift_indices):
            if i==0: plt.axvline(x=rmrs[v], color='k', linestyle='dotted', label='$\Delta\alpha_B$')
            else: plt.axvline(x=rmrs[v], color='k', linestyle='dotted')
    
    plt.plot(rmrsPlotting, avg_flight_angle, 'darkkhaki')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flight Angle [rad]')
    plt.title('Average Flight Angle before Striking the Surface')
    #plt.legend()
    plt.show(block=plot_blocker)
    #if not plot_blocker: plt.savefig('plots/avg_flight_angle.png')
    
    #make histogram of all particles binned into the TIME bins
    histogram_particle_time_AllSegs = np.sum(histogram_particle_time, axis=0)
    nP_time = np.sum(histogram_particle_time_AllSegs)
    print('Time nP:', nP_time)
    
    plt.close()
    plt.bar(bin_edges_time[:-1], histogram_particle_time_AllSegs/nP_time, width=bin_width_time, color='darkviolet', align='edge', edgecolor='k')
    plt.axvline(np.log10(3e-6), 0,1, color='darkkhaki')
    plt.axvline(np.log10(50e-6), 0,1, color='darkkhaki')
    plt.xlabel('Logarithmic Time [log(sec)]')
    plt.ylabel('Counts\n')
    plt.title('Flight Time before Striking the Surface')# \n nP=%.4E'%(nP))
    plt.show(block=plot_blocker)
    
    #make histogram of all particles binned into the ANGLE bins
    histogram_particle_angle_AllSegs = np.sum(histogram_particle_angle, axis=0)
    nP_angle = np.sum(histogram_particle_angle_AllSegs)
    print('Angle nP:', nP_angle)
    
    plt.close()
    plt.bar(bin_edges_angle[:-1], histogram_particle_angle_AllSegs/nP_angle, width=bin_width_angle, color='darkkhaki', align='edge', edgecolor='k')
    plt.axvline(np.log10(2*np.pi), 0,1, color='darkviolet')
    plt.xlabel('Logarithmic Angle [log(rad)]')
    plt.ylabel('Counts\n')
    plt.title('Flight Angle before Striking the Surface')# \n nP=%.4E'%(nP))
    plt.show(block=plot_blocker)
    
    ########################################################
    # Calculate percentage under a peak
    ########################################################
    
    cutoff = -5.5
    peak_edge_time = np.where((np.abs(bin_edges_time - cutoff) == np.min(np.abs(bin_edges_time - cutoff))))[0][0]
    prompt_redep_time = np.sum(histogram_particle_time_AllSegs[:round(peak_edge_time)])/nP_time
    print('Fraction Promptly Redeposited by Flight Time:', prompt_redep_time)
    
    cutoff = np.log10(2*np.pi)
    peak_edge_angle = np.where((np.abs(bin_edges_angle - cutoff) == np.min(np.abs(bin_edges_angle - cutoff))))[0][0]
    prompt_redep_angle = np.sum(histogram_particle_angle_AllSegs[:round(peak_edge_angle)])/nP_angle
    print('Fraction Promptly Redeposited by Flight Angle:', prompt_redep_angle)
        
    return


if __name__ == "__main__":
    #plot_history2D(run_directory+'/output/history.nc')
    plot_surf_nc([1,6], 9, [1,6], '../examples/sasvw-pa-fav/output/perlmutter/production/surface_S1.nc', \
                 '../examples/sasvw-pa-fav/output/perlmutter/production/positions_S1.nc')
    #plot_surf_nc([1,4], 9, [1,5], '../examples/sasvw-pa-fav/output/surface1.nc', \
                 #'../examples/sasvw-pa-fav/output/positions1.nc')
    #plot_surf_nc([5,2], 8, [1,5], setup_directory+"/../output/perlmutter/production/forces24.09.19/surfaces/BFT.nc", \
                 #setup_directory+'/../output/perlmutter/production/forces24.09.19/positions/BFT.nc', norm='')
    #analyze_leakage('perlmutter/history_D3t6.nc')
    #analyze_forces('ExB drift', 'z', rzlim=True, colorbarLimits=[-500,500], dt=1e-9)
    
    #init()
    #plot_gitr_gridspace()
    #plot_particle_source()
    #plot_history2D(setup_directory+"/../output/perlmutter/production/forces24.09.19/histories/gradT.nc",\
    #plot_history2D(setup_directory+"/../output/history1.nc",\
    #plot_history2D(setup_directory+"/../output/history1.nc",\
    #plot_history2D("/pscratch/sd/h/hayes/sasvw-pa-fav-history/output/history.nc",\
                   #bFile=setup_directory+'/../input/bField.nc')
    #spectroscopy(2013859273149157.8,3,specFile='/Users/Alyssa/Desktop/spec.nc')
    #ionization_analysis([0,0], '../examples/sasvw-pa-fav/output/perlmutter/production/','history_IF.nc', 'positions_IF.nc')
    #prompt_redep_hist([2,8,5], '../examples/sasvw-pa-fav/output/perlmutter/production/forces24.09.19/positions/','BEF.nc')
    #particle_diagnostics_hist(4, '/Users/Alyssa/Dev/GITR_processing/examples/sasvw-pa-fav/output/perlmutter/production/particle_histograms_PR1.nc', plot_blocker=True)
    #particle_diagnostics_hist(1e4, run_directory+'/output/particle_histograms.nc', plot_blocker=False)
