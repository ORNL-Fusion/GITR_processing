import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4

################################################
# setting directories and special constants
################################################

setup_directory = '.'
run_directory = '/pscratch/sd/h/hayes/sasvw-vertex-unfav/history'
rmrs_fine_file = setup_directory+'/assets/rmrs_fine.txt'

W_surf_indices = np.arange(16,25)
tile_shift_indices = [2,6]
Bangle_shift_indices = [3,6]
r_sp, z_sp = 1.50230407, 1.23187366 #vertex

sys.path.insert(0, os.path.abspath(setup_directory))
sys.path.insert(0, os.path.abspath('../../../python'))

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
    
    plt.xlim(1.37, 1.52)
    plt.ylim(1.06, 1.23)
    plt.title('Case 4: W Trajectories', fontsize=20)
    #plt.show(block=False)
    plt.savefig('history.svg')
    plt.close()
    return


if __name__ == "__main__":
    plot_history2D(run_directory+"/output/history.nc",\
                   bFile=setup_directory+'/../input/bField.nc')
