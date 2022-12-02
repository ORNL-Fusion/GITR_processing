import numpy as np
import netCDF4
import matplotlib.pyplot as plt

#import wall geometry to plot over
wall_file = '../setup/assets/gitr_rz.txt'
with open(wall_file, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])

#import W surface indices
surf_file = '../setup/assets/surf_ind.txt'
with open(surf_file, 'r') as file:
    surf = file.readlines()

surf = np.array(surf,dtype='int')
Rsurf = np.zeros(len(surf))
Zsurf = np.zeros(len(surf))
for i,v in enumerate(surf):
    Rsurf[i] = R[v]
    Zsurf[i] = Z[v]

area = 0.12*np.sqrt(np.power((Rsurf[:-1]-Rsurf[1:]),2) + np.power((Zsurf[:-1]-Zsurf[1:]),2))
dist = np.sqrt((Rsurf[:-1]-Rsurf[1:])**2 + (Zsurf[:-1]-Zsurf[1:])**2)
rmrs =  np.cumsum(dist)

def plot_history2D():
    history = netCDF4.Dataset("history.nc", "r", format="NETCDF4")

    x = history.variables['x']
    z = history.variables['z']

    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Particle Tracks')

    for i in range(0,len(x)):
        plt.plot(x[i][:],z[i][:], linewidth=0.2)

    plt.xlim(1.43, 1.56)
    plt.ylim(1.07, 1.25)
    plt.savefig('plots/history.pdf')

def plot_surf_nc(ls,fs, pps_per_nP):
    surface = netCDF4.Dataset("surface.nc", "r", format="NETCDF4")
    print(surface.variables['grossErosion'][:])

    grossEro = (surface.variables['grossErosion'][:])*pps_per_nP
    grossDep = (surface.variables['grossDeposition'][:])*pps_per_nP
    netEro = (grossEro-grossDep)
    
    grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)/area
    grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)/area
    netEro = np.average([netEro[:-1], netEro[1:]],axis=0)/area
    
    print('rmrs',len(rmrs))
    print('surf',len(grossEro))

    plt.close()
    plt.plot(rmrs,grossEro,'r', label='Gross Erosion',linewidth=ls)
    plt.plot(rmrs,grossDep,'g', label='Redeposition',linewidth=ls)
    plt.plot(rmrs,netEro,'k', label='Net Erosion',linewidth=ls)
    plt.xlabel('D-Dsep [m]',fontsize=fs)
    plt.ylabel('Flux [#/m2s]',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(loc='upper left', fontsize=fs)
    plt.title('GITR Predicted Erosion and Redeposition Profiles',fontsize=fs)
    plt.savefig('plots/surface.pdf')

def plot_surf_plasma_params(plot_variables):
    profilesFile = '../input/profiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    r_mesh = profiles.variables['r'][:]
    z_mesh = profiles.variables['z'][:]
    
    r_mid = np.average(np.array([Rsurf[:-1],Rsurf[1:]]),axis=0)
    z_mid = np.average(np.array([Zsurf[:-1],Zsurf[1:]]),axis=0)
    r_indices = interpolate(r_mid,r_mesh)
    z_indices = interpolate(z_mid,z_mesh)
    
    #scoot the r_indices 1 cell to the left if the profiles.nc gridcell is too far off
    br = profiles.variables['br'][:][z_indices,r_indices]
    for i in range(len(br)):
        if br[i]==-1: r_indices[i] -= 1
    
    ne = profiles.variables['ne'][:][z_indices, r_indices]
    te = profiles.variables['te'][:][z_indices, r_indices]
    
    if plot_variables==1:    
        ls, fs = 3, 14
        plt.rcParams.update({'font.size':fs})
        
        #electron density along surface
        plt.close()
        plt.plot(rmrs, ne, linewidth=ls)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('ne [m-3]')
        plt.title('Electron Density along the SAS-VW Divertor')
        plt.savefig('plots/surf_ne.png')
    
        #electron temperature along surface
        plt.close()
        plt.plot(rmrs, te, linewidth=ls)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Te [eV]')
        plt.title('Electron Temperature along the SAS-VW Divertor')
        plt.savefig('plots/surf_te.png')
        
        #electron temperature dependence of ionization rate
        plt.close()
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
        plt.plot(kte, invtau0, 'red', linewidth=ls, label='0 → 1')
        plt.plot(kte, invtau1, 'darkorange', linewidth=ls, label='1 → 2')
        plt.plot(kte, invtau2, 'gold', linewidth=ls, label='2 → 3')
        plt.plot(kte, invtau3, 'green', linewidth=ls, label='3 → 4')
        plt.plot(kte, invtau4, 'blue', linewidth=ls, label='4 → 5')
        plt.xlabel('Te [eV]')
        plt.ylabel('Relative Ionization Rate [1/s]')
        plt.title('Te Dependence of Ionization Rate for W')
        plt.legend()
        plt.savefig('plots/IonizRate(Te)')
        
        #ionization rate along surface
        surf_IonizRate0 = ne/(np.exp(x0/te)/np.sqrt(te))
        surf_IonizRate1 = ne/(np.exp(x1/te)/np.sqrt(te))
        surf_IonizRate2 = ne/(np.exp(x2/te)/np.sqrt(te))
        surf_IonizRate3 = ne/(np.exp(x3/te)/np.sqrt(te))
        surf_IonizRate4 = ne/(np.exp(x4/te)/np.sqrt(te))
        plt.close()
        plt.plot(rmrs, surf_IonizRate0, 'red', linewidth=ls, label='0 → 1')
        plt.plot(rmrs, surf_IonizRate1, 'darkorange', linewidth=ls, label='1 → 2')
        plt.plot(rmrs, surf_IonizRate2, 'gold', linewidth=ls, label='2 → 3')
        plt.plot(rmrs, surf_IonizRate3, 'green', linewidth=ls, label='3 → 4')
        plt.plot(rmrs, surf_IonizRate4, 'blue', linewidth=ls, label='4 → 5')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Relative Ionization Rate [1/s]')
        plt.title('Ionization Rate along SAS-VW')
        plt.legend()
        plt.savefig('plots/surf_IonizRate')

        
        

def interpolate(small,big):
    indices = np.zeros(len(small))
    for i in range(len(small)):
        diff = np.min(np.abs(big-small[i]))
        index_possibilities = np.array([np.where(big==small[i]+diff)[0], \
                             np.where(big==small[i]-diff)[0]], dtype=object)
        try: indices[i] = index_possibilities[1]
        except: indices[i] = index_possibilities[0]

    return indices.astype(int)



if __name__ == "__main__":
    #plot_history2D()
    #plot_surf_nc(1,30,111879178639.80714)
    plot_surf_plasma_params(1)
